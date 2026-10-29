//! Pauli-network synthesis: a circuit of Clifford gates and Z rotations
//! rewritten as Pauli rotations, then synthesized again from scratch.
//!
//! Every Z rotation of the circuit, moved through the Clifford gates before
//! it, is a rotation `exp(-iθ/2 P)` about some Pauli string `P`; the circuit
//! is these rotations, in order, followed by one Clifford. The synthesis
//! (after Rustiq, Goubault de Brugière and Martiel, "Faster and shorter
//! synthesis of Hamiltonian simulation circuits", 2024, and pytket's
//! GreedyPauliSimp) builds the circuit again: it applies the two-qubit
//! Clifford gate that most lowers the weights of the rotations that can come
//! next, and emits a rotation as an Rz once it acts on one qubit. The
//! Clifford left at the end is synthesized by decoupling one qubit at a time,
//! as in the greedy synthesis of Bravyi, Shaydulin, Hu and Maslov
//! ("Clifford circuit optimization with templates and symbolic Pauli
//! gates", 2021).
//!
//! Single-qubit gates are most of what such networks cost, so three choices
//! keep them few:
//!
//! - The local Cliffords that set up a CX are not undone after it. The
//!   rotations still to come are read in the new frame, which changes none
//!   of their weights, and the frame is paid for once, by the Clifford left
//!   at the end. Undoing them nearly doubled the single-qubit gates.
//! - Each CX is placed whichever way round needs fewer local Cliffords: X on
//!   `a` controlling Z on `b` is a bare CX(b, a), where CX(a, b) needs a
//!   basis change on both qubits.
//! - Basis changes use √X where a Hadamard would do (Y to Z is one √X; X to
//!   Z is S then √X), as √X is native on IBM-style devices and H costs three
//!   gates there.
//!
//! Commuting rotations may be reordered, and CNOTs shared between them,
//! which sequences of rotations written one after the other cannot do: on
//! the QFT this finds a third fewer CNOTs. It is a candidate the pipeline
//! keeps only when it needs no more two-qubit gates and no more gates in
//! all.

use crate::ir::Instruction;
use crate::rng::Rng;
use std::f64::consts::FRAC_PI_2;

/// Circuits up to this many qubits (bits of the Pauli bitsets).
pub const MAX_QUBITS: usize = 128;
/// Circuits with more rotations than this are not worth the time.
const MAX_ROTATIONS: usize = 400;
/// Nor circuits with fewer two-qubit gates per rotation than this. A Pauli
/// exponential's CX ladders spend 2(w − 1) on a rotation of weight w; QFTs,
/// adders, Clifford+T and QAOA about one, and the network never beat the
/// direct compilation on those.
const CX_PER_ROTATION: usize = 3;

/// Whether the synthesis is worth trying on `instructions`, by their
/// two-qubit gates per non-Clifford Z rotation.
pub fn promising(instructions: &[Instruction]) -> bool {
    let mut rotations = 0;
    let mut two_qubit = 0;
    for inst in instructions {
        match (inst.name.as_str(), inst.params.as_slice()) {
            (_, _) if inst.qubits.len() == 2 => two_qubit += 1,
            ("t" | "tdg", []) => rotations += 1,
            ("rz" | "p" | "u1" | "phase", &[theta]) => {
                let quarter = theta / FRAC_PI_2;
                if (quarter - quarter.round()).abs() > 1e-12 {
                    rotations += 1;
                }
            }
            _ => {}
        }
    }
    rotations > 0 && rotations <= MAX_ROTATIONS && two_qubit >= CX_PER_ROTATION * rotations
}

/// `i^phase · Π_q X_q^{x_q} Z_q^{z_q}`, each qubit's X before its Z.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Pauli {
    x: u128,
    z: u128,
    phase: u8,
}

impl Pauli {
    const IDENTITY: Pauli = Pauli {
        x: 0,
        z: 0,
        phase: 0,
    };

    fn x(q: usize) -> Pauli {
        Pauli {
            x: 1 << q,
            ..Pauli::IDENTITY
        }
    }

    fn z(q: usize) -> Pauli {
        Pauli {
            z: 1 << q,
            ..Pauli::IDENTITY
        }
    }

    fn weight(&self) -> u32 {
        (self.x | self.z).count_ones()
    }

    fn support(&self) -> u128 {
        self.x | self.z
    }

    /// The letter on qubit q: 0 = I, 1 = X, 2 = Z, 3 = Y (as X and Z bits).
    fn letter(&self, q: usize) -> u8 {
        (((self.x >> q) & 1) as u8) | ((((self.z >> q) & 1) as u8) << 1)
    }

    fn commutes(&self, other: &Pauli) -> bool {
        ((self.x & other.z) ^ (self.z & other.x))
            .count_ones()
            .is_multiple_of(2)
    }

    /// `self · other`.
    fn times(&self, other: &Pauli) -> Pauli {
        // Moving other's X past self's Z on the same qubit gives a −1.
        let swaps = (self.z & other.x).count_ones() as u8;
        Pauli {
            x: self.x ^ other.x,
            z: self.z ^ other.z,
            phase: (self.phase + other.phase + 2 * swaps) % 4,
        }
    }

    /// The real sign of a Hermitian Pauli: its phase once each Y (= i X Z) is
    /// counted, ±1.
    fn sign(&self) -> f64 {
        let ys = (self.x & self.z).count_ones() as u8;
        // i^phase X^x Z^z = i^(phase − ys) · Π (Y for each XZ pair, …)
        if (self.phase + 4 - ys % 4).is_multiple_of(4) {
            1.0
        } else {
            -1.0
        }
    }
}

/// A Clifford gate the synthesis uses: on one qubit, or a CX.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Gate {
    H(usize),
    S(usize),
    Sdg(usize),
    /// √X, which is H S H.
    Sx(usize),
    Cx(usize, usize),
}

impl Gate {
    /// `g P g†`.
    fn conjugate(&self, p: &mut Pauli) {
        match *self {
            Gate::H(q) => {
                let (x, z) = ((p.x >> q) & 1, (p.z >> q) & 1);
                p.phase = (p.phase + 2 * (x & z) as u8) % 4;
                p.x = (p.x & !(1 << q)) | (z << q);
                p.z = (p.z & !(1 << q)) | (x << q);
            }
            Gate::S(q) => {
                let x = (p.x >> q) & 1;
                p.phase = (p.phase + x as u8) % 4;
                p.z ^= x << q;
            }
            Gate::Sdg(q) => {
                let x = (p.x >> q) & 1;
                p.phase = (p.phase + 3 * x as u8) % 4;
                p.z ^= x << q;
            }
            Gate::Sx(q) => {
                for g in [Gate::H(q), Gate::S(q), Gate::H(q)] {
                    g.conjugate(p);
                }
            }
            Gate::Cx(c, t) => {
                p.x ^= ((p.x >> c) & 1) << t;
                p.z ^= ((p.z >> t) & 1) << c;
            }
        }
    }

    /// The inverse, for the gates [`read`] gives.
    fn inverse(&self) -> Gate {
        match *self {
            Gate::Sx(_) => unreachable!("read gives no √X"),
            Gate::S(q) => Gate::Sdg(q),
            Gate::Sdg(q) => Gate::S(q),
            g => g,
        }
    }

    fn instruction(&self) -> Instruction {
        match *self {
            Gate::H(q) => Instruction::gate("h", vec![q], vec![]),
            Gate::S(q) => Instruction::gate("s", vec![q], vec![]),
            Gate::Sdg(q) => Instruction::gate("sdg", vec![q], vec![]),
            Gate::Sx(q) => Instruction::gate("sx", vec![q], vec![]),
            Gate::Cx(c, t) => Instruction::gate("cx", vec![c, t], vec![]),
        }
    }
}

/// The circuit's gates as [`Gate`]s and Z rotations (`None` for a gate this
/// module does not handle).
fn read(inst: &Instruction) -> Option<Vec<Step>> {
    let q = |i: usize| inst.qubits[i];
    let gates = |gs: &[Gate]| Some(gs.iter().map(|&g| Step::Clifford(g)).collect());
    if !inst.is_gate() {
        return None;
    }
    match (inst.name.as_str(), inst.params.as_slice()) {
        ("h", []) => gates(&[Gate::H(q(0))]),
        ("s", []) => gates(&[Gate::S(q(0))]),
        ("sdg", []) => gates(&[Gate::Sdg(q(0))]),
        ("z", []) => gates(&[Gate::S(q(0)), Gate::S(q(0))]),
        ("x", []) => gates(&[Gate::H(q(0)), Gate::S(q(0)), Gate::S(q(0)), Gate::H(q(0))]),
        ("y", []) => gates(&[
            Gate::S(q(0)),
            Gate::S(q(0)),
            Gate::H(q(0)),
            Gate::S(q(0)),
            Gate::S(q(0)),
            Gate::H(q(0)),
        ]),
        ("sx", []) => gates(&[Gate::H(q(0)), Gate::S(q(0)), Gate::H(q(0))]),
        ("sxdg", []) => gates(&[Gate::H(q(0)), Gate::Sdg(q(0)), Gate::H(q(0))]),
        ("cx" | "cnot", []) => gates(&[Gate::Cx(q(0), q(1))]),
        ("cz", []) => gates(&[Gate::H(q(1)), Gate::Cx(q(0), q(1)), Gate::H(q(1))]),
        ("swap", []) => gates(&[
            Gate::Cx(q(0), q(1)),
            Gate::Cx(q(1), q(0)),
            Gate::Cx(q(0), q(1)),
        ]),
        ("rz" | "p" | "u1" | "phase", &[theta]) => Some(vec![Step::Rotation(q(0), theta)]),
        ("t", []) => Some(vec![Step::Rotation(q(0), std::f64::consts::FRAC_PI_4)]),
        ("tdg", []) => Some(vec![Step::Rotation(q(0), -std::f64::consts::FRAC_PI_4)]),
        ("id", []) => Some(vec![]),
        _ => None,
    }
}

/// A rotation `exp(-iθ/2 P)`, as `(P, θ)`.
type Rotation = (Pauli, f64);
/// A rotation and the index of the instruction it comes from.
type Located = (Pauli, f64, usize);

enum Step {
    Clifford(Gate),
    Rotation(usize, f64),
}

/// The circuit as Pauli rotations, in order, and the inverse of the Clifford
/// that follows them, as the images `C† X_q C` and `C† Z_q C`; None if it
/// has a gate other than Clifford gates and Z rotations, or too many qubits.
fn extract(instructions: &[Instruction], n: usize) -> Option<(Vec<Rotation>, Vec<Pauli>)> {
    let (rotations, m) = extract_indexed(instructions, n)?;
    Some((
        rotations
            .into_iter()
            .map(|(p, theta, _)| (p, theta))
            .collect(),
        m,
    ))
}

/// [`extract`], with the index of the instruction each rotation comes from.
fn extract_indexed(instructions: &[Instruction], n: usize) -> Option<(Vec<Located>, Vec<Pauli>)> {
    if n > MAX_QUBITS {
        return None;
    }
    // m[2q] = C† X_q C and m[2q + 1] = C† Z_q C, for the Clifford C so far.
    let mut m: Vec<Pauli> = (0..n).flat_map(|q| [Pauli::x(q), Pauli::z(q)]).collect();
    let image = |m: &[Pauli], p: &Pauli| -> Pauli {
        let mut out = Pauli {
            phase: p.phase,
            ..Pauli::IDENTITY
        };
        for q in 0..n {
            if (p.x >> q) & 1 == 1 {
                out = out.times(&m[2 * q]);
            }
            if (p.z >> q) & 1 == 1 {
                out = out.times(&m[2 * q + 1]);
            }
        }
        out
    };
    let mut rotations = Vec::new();
    for (index, inst) in instructions.iter().enumerate() {
        for step in read(inst)? {
            match step {
                Step::Rotation(q, theta) => {
                    // Multiples of π/2 are Clifford: S, Z or S†.
                    let quarter = theta / FRAC_PI_2;
                    if (quarter - quarter.round()).abs() < 1e-12 {
                        let k = quarter.round().rem_euclid(4.0) as usize;
                        for _ in 0..k {
                            apply_to_map(&mut m, Gate::S(q), &image);
                        }
                    } else {
                        rotations.push((m[2 * q + 1], theta, index));
                    }
                }
                Step::Clifford(g) => apply_to_map(&mut m, g, &image),
            }
        }
    }
    Some((rotations, m))
}

/// Merge Z rotations that act about the same Pauli string once moved through
/// the Clifford gates before them, when every rotation between them
/// commutes with it: the later one takes the sum and the earlier goes.
///
/// This is phase folding through any Clifford gate, Hadamards included,
/// after Zhang and Chen ("Optimizing T gates in Clifford+T circuit as π/4
/// rotations around Paulis", 2019). Circuits with other gates (or more than
/// [`MAX_QUBITS`] qubits) are left alone.
pub fn merge_rotations(instructions: Vec<Instruction>) -> Vec<Instruction> {
    let n = instructions
        .iter()
        .flat_map(|i| i.qubits.iter().map(|&q| q + 1))
        .max()
        .unwrap_or(0);
    let Some((rotations, _)) = extract_indexed(&instructions, n) else {
        return instructions;
    };
    // The rotations still standing, by position in `rotations`, and each
    // one's angle about its Pauli string with the sign taken out.
    let mut alive: Vec<usize> = Vec::new();
    let mut angle: Vec<f64> = rotations.iter().map(|(p, t, _)| p.sign() * t).collect();
    let mut gone = vec![false; rotations.len()];
    let mut changed = false;
    for k in 0..rotations.len() {
        let pk = rotations[k].0;
        let mut i = alive.len();
        let mut looked = 0;
        while i > 0 && looked < MERGE_LOOKBACK {
            i -= 1;
            looked += 1;
            let j = alive[i];
            let pj = rotations[j].0;
            if (pj.x, pj.z) == (pk.x, pk.z) {
                angle[k] += angle[j];
                gone[j] = true;
                alive.remove(i);
                changed = true;
            } else if !pj.commutes(&pk) {
                break;
            }
        }
        alive.push(k);
    }
    if !changed {
        return instructions;
    }
    // Each instruction's new angle (or removal), from its rotation's.
    let mut removed = vec![false; instructions.len()];
    let mut new_angle: Vec<Option<f64>> = vec![None; instructions.len()];
    for (k, &(p, theta, index)) in rotations.iter().enumerate() {
        let value = angle[k] * p.sign();
        let turns = value.rem_euclid(std::f64::consts::TAU);
        if gone[k] || turns < 1e-12 || std::f64::consts::TAU - turns < 1e-12 {
            removed[index] = true;
        } else if (value - theta).abs() > 1e-12 {
            new_angle[index] = Some(value);
        }
    }
    instructions
        .into_iter()
        .enumerate()
        .filter(|(i, _)| !removed[*i])
        .map(|(i, inst)| match new_angle[i] {
            Some(theta) => Instruction::gate("rz", inst.qubits, vec![theta]),
            None => inst,
        })
        .collect()
}

/// How many rotations [`merge_rotations`] looks back past, at most.
const MERGE_LOOKBACK: usize = 256;

/// `m ← m ∘ conj(g†)`: appending `g` to the Clifford C makes
/// `(gC)† P (gC) = C† (g† P g) C`.
fn apply_to_map(m: &mut [Pauli], g: Gate, image: &dyn Fn(&[Pauli], &Pauli) -> Pauli) {
    let support: &[usize] = match g {
        Gate::H(q) | Gate::S(q) | Gate::Sdg(q) | Gate::Sx(q) => &[q, q],
        Gate::Cx(c, t) => &[c, t],
    };
    let inverse = g.inverse();
    let mut updates = Vec::with_capacity(4);
    for &q in support {
        for (row, generator) in [(2 * q, Pauli::x(q)), (2 * q + 1, Pauli::z(q))] {
            let mut p = generator;
            inverse.conjugate(&mut p);
            updates.push((row, image(m, &p)));
        }
    }
    for (row, p) in updates {
        m[row] = p;
    }
}

/// How far ahead of the next rotation the synthesis looks.
const WINDOW: usize = 48;
/// How much each later commuting layer counts, relative to the one before.
const DISCOUNT: f64 = 0.3;
/// How much the Clifford left counts, per qubit its images act on beyond one.
const TABLEAU_WEIGHT: f64 = 0.1;
/// What each single-qubit Clifford of a candidate gate costs, next to the
/// weights it removes: enough to break ties.
const LOCAL_PENALTY: f64 = 0.01;
/// Pairs of qubits tried as places for the next gate, at most.
const MAX_PAIRS: usize = 64;
/// With seed k, at most k times this much random noise is added to each
/// candidate gate's score. Well under a weight, it only reorders gates that
/// do about as well, which is where the greedy synthesis is blind; the seeds
/// spread from gentle to bolder reorderings.
const NOISE: f64 = 0.05;

/// Local Cliffords on `q` taking `letter` (1 = X, 2 = Z, 3 = Y) to Z.
fn to_z(q: usize, letter: u8) -> Vec<Gate> {
    match letter {
        1 => vec![Gate::S(q), Gate::Sx(q)],
        3 => vec![Gate::Sx(q)],
        _ => vec![],
    }
}

/// Local Cliffords on `q` taking `letter` to X.
fn to_x(q: usize, letter: u8) -> Vec<Gate> {
    match letter {
        2 => vec![Gate::Sx(q), Gate::S(q)],
        3 => vec![Gate::Sdg(q)],
        _ => vec![],
    }
}

/// The circuit synthesized again as a Pauli network, and where each qubit's
/// state ends (`wires[q]`; always `q` unless `permute`), or None if it has
/// gates other than Clifford gates and Z rotations (or more than
/// [`MAX_QUBITS`] qubits).
pub fn synthesize(
    instructions: &[Instruction],
    permute: bool,
) -> Option<(Vec<Instruction>, Vec<usize>)> {
    synthesize_seeded(instructions, permute, 0)
}

/// [`synthesize`], breaking ties between candidate gates at random when
/// `seed` is not 0. The greedy choice of each CX looks only a few layers
/// ahead, and between gates that do equally well there it can lose a tenth of
/// the CNOTs; a few seeds, run side by side, find the better ones.
pub fn synthesize_seeded(
    instructions: &[Instruction],
    permute: bool,
    seed: u64,
) -> Option<(Vec<Instruction>, Vec<usize>)> {
    let n = instructions
        .iter()
        .flat_map(|i| i.qubits.iter().map(|&q| q + 1))
        .max()
        .unwrap_or(0);
    let (rotations, rows) = extract(instructions, n)?;

    let mut state = State {
        pending: rotations,
        rows,
        out: Vec::new(),
        n,
        rng: (seed != 0).then(|| Rng::new(seed)),
        noise: NOISE * seed as f64,
    };
    state.rotations();
    let wires = state.clifford(permute);
    Some((state.out, wires))
}

struct State {
    /// The rotations still to emit, in the current frame.
    pending: Vec<Rotation>,
    /// The Clifford still to emit, inverted: `rows[2q] = E X_q E†`, …
    rows: Vec<Pauli>,
    out: Vec<Instruction>,
    n: usize,
    /// Breaks ties between candidate gates, when seeded, adding up to
    /// `noise` to each one's score.
    rng: Option<Rng>,
    noise: f64,
}

impl State {
    /// Emit `g` first: every rotation still to come, and the rest, move
    /// through it.
    fn emit(&mut self, g: Gate) {
        for (p, _) in self.pending.iter_mut() {
            g.conjugate(p);
        }
        for p in self.rows.iter_mut() {
            g.conjugate(p);
        }
        self.out.push(g.instruction());
    }

    /// Turn the letter on `q` into Z.
    fn turn_to_z(&mut self, letter: u8, q: usize) {
        for g in to_z(q, letter) {
            self.emit(g);
        }
    }

    /// The commuting layer of each rotation among the first [`WINDOW`]: one
    /// more than the latest earlier rotation it does not commute with.
    fn layers(&self) -> Vec<usize> {
        let window = self.pending.len().min(WINDOW);
        let mut layer = vec![0usize; window];
        for i in 0..window {
            for j in 0..i {
                if !self.pending[i].0.commutes(&self.pending[j].0) {
                    layer[i] = layer[i].max(layer[j] + 1);
                }
            }
        }
        layer
    }

    fn rotations(&mut self) {
        let mut stuck = 0;
        while !self.pending.is_empty() {
            let layers = self.layers();
            // Emit a rotation of the first layer that acts on one qubit.
            let ready =
                (0..layers.len()).find(|&i| layers[i] == 0 && self.pending[i].0.weight() <= 1);
            if let Some(i) = ready {
                let (p, theta) = self.pending[i];
                if p.weight() == 1 {
                    let q = p.support().trailing_zeros() as usize;
                    self.turn_to_z(p.letter(q), q);
                    let (p, _) = self.pending[i];
                    let angle = p.sign() * theta;
                    self.out.push(Instruction::gate("rz", vec![q], vec![angle]));
                }
                self.pending.remove(i);
                stuck = 0;
                continue;
            }
            let gates = match stuck < 4 * self.n {
                true => self.best_gate(&layers),
                false => None,
            }
            .unwrap_or_else(|| {
                let first = (0..layers.len())
                    .find(|&i| layers[i] == 0)
                    .expect("a first layer");
                self.reduce_first(first)
            });
            for g in gates {
                self.emit(g);
            }
            stuck += 1;
        }
    }

    /// The Paulis whose weights a gate on qubits `a` and `b` can change: the
    /// rotations in the window, discounted by layer, and the Clifford's
    /// images, weighted by [`TABLEAU_WEIGHT`]. Each as (how much it counts,
    /// its weight, its letters on `a` and `b`).
    fn touching(&self, layers: &[usize], a: usize, b: usize) -> Vec<(f64, u32, usize)> {
        let mask = (1u128 << a) | (1u128 << b);
        let entry = |count: f64, p: &Pauli| {
            (
                count,
                p.weight(),
                (p.letter(a) | (p.letter(b) << 2)) as usize,
            )
        };
        let rotations = layers
            .iter()
            .enumerate()
            .map(|(i, &layer)| (DISCOUNT.powi(layer as i32), &self.pending[i].0));
        let images = self.rows.iter().map(|p| (TABLEAU_WEIGHT, p));
        rotations
            .chain(images)
            .filter(|(_, p)| p.support() & mask != 0)
            .map(|(count, p)| entry(count, p))
            .collect()
    }

    /// How much `gates`, on qubits `a` and `b`, change the cost of the
    /// rotations in the window and of the Clifford left: the weights beyond
    /// one of the Paulis [`State::touching`] gives. Only their letters on
    /// `a` and `b` change, and only as those letters say, so the gates are
    /// applied once to each of the 16 pairs of letters rather than to every
    /// Pauli.
    fn change(gates: &[Gate], touching: &[(f64, u32, usize)], a: usize, b: usize) -> f64 {
        let mut letters_after = [0u32; 16];
        for (pair, after) in letters_after.iter_mut().enumerate() {
            let mut p = Pauli::IDENTITY;
            for (q, letter) in [(a, pair & 3), (b, pair >> 2)] {
                p.x |= ((letter & 1) as u128) << q;
                p.z |= ((letter >> 1) as u128) << q;
            }
            gates.iter().for_each(|g| g.conjugate(&mut p));
            *after = p.weight();
        }
        let letters_before = |pair: usize| u32::from(pair & 3 != 0) + u32::from(pair >> 2 != 0);
        touching
            .iter()
            .map(|&(count, weight, pair)| {
                let after = weight - letters_before(pair) + letters_after[pair];
                count * (f64::from(after.max(1)) - f64::from(weight.max(1)))
            })
            .sum()
    }

    /// The controlled-Pauli gate, σ on `a` controlling τ on `b`, as local
    /// Cliffords then a CX (which is Z controlling X), whichever way round
    /// needs fewer local Cliffords. The local Cliffords are not undone after
    /// the CX (see the module's documentation), so different ways round leave
    /// different frames, but the same weights.
    fn controlled_pauli(sigma: u8, tau: u8, a: usize, b: usize) -> Vec<Gate> {
        let [forward, backward] = State::both_ways(sigma, tau, a, b);
        if backward.len() < forward.len() {
            backward
        } else {
            forward
        }
    }

    /// σ on `a` controlling τ on `b` with a CX(a, b), and as τ on `b`
    /// controlling σ on `a`, with a CX(b, a).
    fn both_ways(sigma: u8, tau: u8, a: usize, b: usize) -> [Vec<Gate>; 2] {
        [
            [to_z(a, sigma), to_x(b, tau), vec![Gate::Cx(a, b)]].concat(),
            [to_x(a, sigma), to_z(b, tau), vec![Gate::Cx(b, a)]].concat(),
        ]
    }

    /// The controlled-Pauli gate on two qubits of a first-layer rotation that
    /// lowers [`State::cost`] most, if any lowers it.
    fn best_gate(&mut self, layers: &[usize]) -> Option<Vec<Gate>> {
        let mut pairs = Vec::new();
        for i in (0..layers.len()).filter(|&i| layers[i] == 0) {
            let p = self.pending[i].0;
            let qubits: Vec<usize> = (0..self.n).filter(|&q| p.letter(q) != 0).collect();
            for (k, &a) in qubits.iter().enumerate() {
                for &b in &qubits[k + 1..] {
                    if !pairs.contains(&(a, b)) {
                        pairs.push((a, b));
                    }
                }
            }
            if pairs.len() >= MAX_PAIRS {
                break;
            }
        }
        // (score, change in cost, gates) of the best so far.
        let mut best: Option<(f64, f64, Vec<Gate>)> = None;
        for &(a, b) in &pairs {
            let touching = self.touching(layers, a, b);
            for sigma in 1..=3 {
                for tau in 1..=3 {
                    let gates = State::controlled_pauli(sigma, tau, a, b);
                    let change = State::change(&gates, &touching, a, b);
                    // Of gates that do as well, those with fewer local
                    // Cliffords before their CX.
                    let mut score = change + LOCAL_PENALTY * (gates.len() - 1) as f64;
                    if let Some(rng) = self.rng.as_mut() {
                        score += self.noise * rng.uniform();
                    }
                    if best.as_ref().is_none_or(|(s, _, _)| score < *s - 1e-9) {
                        best = Some((score, change, gates));
                    }
                }
            }
        }
        best.filter(|(_, change, _)| *change < -1e-9)
            .map(|(_, _, g)| g)
    }

    /// Lower the weight of rotation `i` by one, whatever it does to the others.
    fn reduce_first(&self, i: usize) -> Vec<Gate> {
        let p = self.pending[i].0;
        let support = p.support();
        let c = support.trailing_zeros() as usize;
        let t = (support & !(1 << c)).trailing_zeros() as usize;
        State::reduction(&p, c, t)
    }

    /// The gates (local Cliffords, then a CX) that turn the letters on
    /// `(c, t)` into `(I, Z)`.
    fn reduction(p: &Pauli, c: usize, t: usize) -> Vec<Gate> {
        let mut gates = Vec::new();
        for q in [c, t] {
            gates.extend(to_z(q, p.letter(q)));
        }
        gates.push(Gate::Cx(c, t));
        gates
    }

    /// Emit the Clifford left, so that `rows` become the identity, or with
    /// `permute` a permutation of the qubits: one row at a time, the cheapest
    /// first, its pair of images decoupled from the qubits still free and
    /// left on one of them; then the signs. Returns where each row's pair
    /// ends.
    fn clifford(&mut self, permute: bool) -> Vec<usize> {
        let mut rows: Vec<usize> = (0..self.n).collect();
        let mut free: Vec<usize> = (0..self.n).collect();
        let mut wires = vec![0; self.n];
        while !rows.is_empty() {
            let (k, _) = rows
                .iter()
                .enumerate()
                .map(|(k, &q)| (k, self.decoupling_cost(q, &free)))
                .min_by_key(|&(_, c)| c)
                .expect("a row left");
            let q = rows.remove(k);
            let target = self.decouple(q, &free, permute);
            free.retain(|&p| p != target);
            wires[q] = target;
        }
        for q in 0..self.n {
            // rows[2q] = ±X_w, rows[2q+1] = ±Z_w now: fix the signs.
            let w = wires[q];
            if self.rows[2 * q].sign() < 0.0 {
                // Z X Z = −X.
                self.emit(Gate::S(w));
                self.emit(Gate::S(w));
            }
            if self.rows[2 * q + 1].sign() < 0.0 {
                // X Z X = −Z: X = H Z H = H S S H.
                for g in [Gate::H(w), Gate::S(w), Gate::S(w), Gate::H(w)] {
                    self.emit(g);
                }
            }
        }
        wires
    }

    /// Roughly the CXs decoupling row `q` costs: the free qubits its pair of
    /// images acts on, beyond one.
    fn decoupling_cost(&self, q: usize, free: &[usize]) -> u32 {
        let (a, b) = (self.rows[2 * q], self.rows[2 * q + 1]);
        free.iter()
            .map(|&k| match (a.letter(k) != 0, b.letter(k) != 0) {
                (false, false) => 0,
                (true, true) if a.letter(k) == b.letter(k) => 2,
                _ => 1,
            })
            .sum::<u32>()
            .saturating_sub(1)
    }

    /// Make row `q`'s pair `±X_t` and `±Z_t` for a free qubit `t` (`q` itself
    /// unless `permute`), using only the free qubits; returns `t`.
    fn decouple(&mut self, q: usize, free: &[usize], permute: bool) -> usize {
        // The pair anticommutes, so an odd number of qubits carry
        // anticommuting letters; gather them into one, the target.
        let anticommuting = |s: &State, k: usize| {
            let (x, z) = (s.rows[2 * q].letter(k), s.rows[2 * q + 1].letter(k));
            x != 0 && z != 0 && x != z
        };
        let mut odd: Vec<usize> = free
            .iter()
            .copied()
            .filter(|&k| anticommuting(self, k))
            .collect();
        // Prefer q itself as where the pair ends.
        odd.sort_by_key(|&k| k != q);
        // Local Cliffords: anticommuting letters to (X, Z).
        for &k in &odd {
            self.turn_to_xz(q, k);
        }
        // Pairs of the others: CX(k, l) turns (X, Z)(X, Z) into (X, I)(I, Z).
        let target = odd[0];
        let mut rest = odd[1..].to_vec();
        while rest.len() >= 2 {
            let (k, l) = (rest.pop().unwrap(), rest.pop().unwrap());
            self.emit(Gate::Cx(k, l));
        }
        debug_assert!(rest.is_empty());
        // Every other free qubit now carries commuting letters: clear them
        // onto `target`, which carries (X, Z).
        for &k in free.iter().filter(|&&k| k != target) {
            self.clear(q, k, target);
        }
        if target != q && !permute {
            // Move the pair from `target` to `q`, whose letters are I now.
            self.emit(Gate::Cx(target, q));
            self.emit(Gate::Cx(q, target));
            self.emit(Gate::Cx(target, q));
            return q;
        }
        target
    }

    /// Local Cliffords on `k` taking the letters of q's pair on `k`, which
    /// anticommute, to (X, Z).
    fn turn_to_xz(&mut self, q: usize, k: usize) {
        for _ in 0..6 {
            let (x, z) = (self.rows[2 * q].letter(k), self.rows[2 * q + 1].letter(k));
            match (x, z) {
                (1, 2) => return,
                (2, _) | (_, 1) => self.emit(Gate::H(k)),
                _ => self.emit(Gate::S(k)),
            }
        }
        unreachable!("anticommuting letters reach (X, Z)");
    }

    /// Remove qubit `k` from q's pair, whose letters on `target` are (X, Z) and
    /// on `k` commute.
    fn clear(&mut self, q: usize, k: usize, target: usize) {
        let letters = |s: &State| (s.rows[2 * q].letter(k), s.rows[2 * q + 1].letter(k));
        // The same letter on both: Z on both after a local Clifford; then a
        // CZ(target, k) clears X's partner, and CX(k, target) Z's.
        match letters(self) {
            (0, 0) => return,
            (a, b) if a == b => self.turn_to_z(a, k),
            (a, 0) => self.turn_to_z(a, k),
            (0, b) => self.turn_to_z(b, k),
            _ => unreachable!("commuting letters"),
        }
        let (a, b) = letters(self);
        if a == 2 {
            // CZ(target, k) = H(k) CX(target, k) H(k): X_target → X_target Z_k.
            self.emit(Gate::H(k));
            self.emit(Gate::Cx(target, k));
            self.emit(Gate::H(k));
        }
        if b == 2 {
            // CX(k, target): Z_target → Z_k Z_target.
            self.emit(Gate::Cx(k, target));
        }
        debug_assert_eq!(letters(self), (0, 0));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testing::{Rng, circuit_unitary, dense_equal_up_to_phase};

    fn random_circuit(rng: &mut Rng, n: usize, length: usize) -> Vec<Instruction> {
        let names = ["h", "s", "sdg", "x", "sx", "cx", "cz", "rz", "t", "swap"];
        (0..length)
            .map(|_| {
                let name = names[rng.below(names.len())];
                let a = rng.below(n);
                let b = (a + 1 + rng.below(n - 1)) % n;
                match name {
                    "cx" | "cz" | "swap" => Instruction::gate(name, vec![a, b], vec![]),
                    "rz" => Instruction::gate(name, vec![a], vec![rng.range(-3.0, 3.0)]),
                    _ => Instruction::gate(name, vec![a], vec![]),
                }
            })
            .collect()
    }

    #[test]
    fn pauli_products_and_conjugation() {
        // X·Z = −i Y: with Y = i X Z, X Z is i^0 X Z, and Z X = −X Z.
        let (x, z) = (Pauli::x(0), Pauli::z(0));
        assert_eq!(z.times(&x).phase, 2);
        assert!(!x.commutes(&z));
        // H X H = Z; S X S† = Y; CX X_c CX = X_c X_t.
        let mut p = x;
        Gate::H(0).conjugate(&mut p);
        assert_eq!(p, z);
        let mut p = x;
        Gate::S(0).conjugate(&mut p);
        assert_eq!((p.x, p.z, p.sign()), (1, 1, 1.0));
        let mut p = x;
        Gate::Cx(0, 1).conjugate(&mut p);
        assert_eq!((p.x, p.z), (0b11, 0));
    }

    #[test]
    fn random_circuits_keep_their_unitary() {
        let mut rng = Rng::new(8);
        for n in [2, 3, 4] {
            for k in 0..60 {
                let circuit = random_circuit(&mut rng, n, 25);
                for permute in [false, true] {
                    // Seeded or not.
                    let seed = (k % 3) as u64;
                    let (mut out, wires) = synthesize_seeded(&circuit, permute, seed)
                        .expect("Clifford gates and rotations");
                    assert!(
                        out.iter()
                            .all(|i| ["h", "s", "sdg", "sx", "cx", "rz"].contains(&i.name.as_str()))
                    );
                    assert!(permute || wires.iter().enumerate().all(|(q, &w)| q == w));
                    // Move every state back to its own wire.
                    let mut at = wires.clone();
                    for q in 0..n {
                        let w = at[q];
                        if w != q {
                            out.push(Instruction::gate("swap", vec![w, q], vec![]));
                            let p = at.iter().position(|&x| x == q).expect("a permutation");
                            at[p] = w;
                            at[q] = q;
                        }
                    }
                    assert!(
                        dense_equal_up_to_phase(
                            &circuit_unitary(&out, n),
                            &circuit_unitary(&circuit, n),
                            1e-8
                        ),
                        "{circuit:?}"
                    );
                }
            }
        }
    }

    #[test]
    fn controlled_paulis_either_way_change_the_same_weights() {
        let letters = |p: &Pauli| (p.letter(0) != 0, p.letter(1) != 0);
        for sigma in 1..=3 {
            for tau in 1..=3 {
                let [forward, backward] = State::both_ways(sigma, tau, 0, 1);
                // Each Pauli on two qubits keeps the same support either way.
                for pair in 0..16u8 {
                    let p = Pauli {
                        x: u128::from(pair & 1) | (u128::from((pair >> 2) & 1) << 1),
                        z: u128::from((pair >> 1) & 1) | (u128::from(pair >> 3) << 1),
                        phase: 0,
                    };
                    let (mut f, mut b) = (p, p);
                    forward.iter().for_each(|g| g.conjugate(&mut f));
                    backward.iter().for_each(|g| g.conjugate(&mut b));
                    assert_eq!(letters(&f), letters(&b), "σ = {sigma}, τ = {tau}");
                }
                let chosen = State::controlled_pauli(sigma, tau, 0, 1);
                assert_eq!(chosen.len(), forward.len().min(backward.len()));
            }
        }
        // X controlling Z is a CX the other way round, with nothing around it.
        assert_eq!(State::controlled_pauli(1, 2, 0, 1), vec![Gate::Cx(1, 0)]);
        assert_eq!(State::controlled_pauli(2, 1, 0, 1), vec![Gate::Cx(0, 1)]);
    }

    #[test]
    fn rotations_merge_through_hadamards() {
        // T, then H CX H on the other qubit (which commutes with Z on the
        // first), then T again: one rotation of π/2, a Clifford, remains.
        let g = |name: &str, q: Vec<usize>| Instruction::gate(name, q, vec![]);
        let circuit = vec![
            g("t", vec![0]),
            g("h", vec![1]),
            g("cx", vec![0, 1]),
            g("h", vec![1]),
            g("t", vec![0]),
        ];
        let out = merge_rotations(circuit.clone());
        assert_eq!(
            out.iter()
                .filter(|i| i.name == "t" || i.name == "rz")
                .count(),
            1
        );
        assert!(dense_equal_up_to_phase(
            &circuit_unitary(&out, 2),
            &circuit_unitary(&circuit, 2),
            1e-9
        ));
        let mut rng = Rng::new(12);
        let mut removed = 0;
        for _ in 0..200 {
            let n = 3;
            let circuit = random_circuit(&mut rng, n, 30);
            let out = merge_rotations(circuit.clone());
            removed += circuit.len() - out.len();
            assert!(dense_equal_up_to_phase(
                &circuit_unitary(&out, n),
                &circuit_unitary(&circuit, n),
                1e-8
            ));
        }
        assert!(removed > 50, "only {removed} rotations merged");
    }
}
