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
//! gates", 2021); on a device, also with every CX on a coupled pair, each
//! qubit's pair of images gathered along a Steiner tree of the chip, as in
//! Winderl et al. ("Architecture-Aware Synthesis of Stabilizer Circuits
//! from Clifford Tableaus", arXiv:2309.08972).
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
use crate::routing::CouplingGraph;
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
                // H, S, H in one step: X stays, Z becomes i³ X Z (= −Y)
                // and X Z becomes i³ Z (so Y = i X Z becomes Z): the X bit
                // takes the Z bit, and the phase gains three Z bits. `sx_is_h_s_h` checks it against
                // the three gates.
                let z = (p.z >> q) & 1;
                p.phase = (p.phase + 3 * z as u8) % 4;
                p.x ^= z << q;
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
            Gate::H(q) => Instruction::gate("h", [q], []),
            Gate::S(q) => Instruction::gate("s", [q], []),
            Gate::Sdg(q) => Instruction::gate("sdg", [q], []),
            Gate::Sx(q) => Instruction::gate("sx", [q], []),
            Gate::Cx(c, t) => Instruction::gate("cx", [c, t], []),
        }
    }
}

/// The circuit's gates as [`Gate`]s and Z rotations (`None` for a gate this
/// module does not handle).
fn read(inst: &Instruction) -> Option<Step> {
    let q = |i: usize| inst.qubits[i];
    // At most six Cliffords (Y); the rest of the array is never read.
    let gates = |gs: &[Gate]| {
        let mut out = [Gate::H(0); 6];
        out[..gs.len()].copy_from_slice(gs);
        Some(Step::Cliffords(out, gs.len()))
    };
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
        ("rz" | "p" | "u1" | "phase", &[theta]) => Some(Step::Rotation(q(0), theta)),
        ("t", []) => Some(Step::Rotation(q(0), std::f64::consts::FRAC_PI_4)),
        ("tdg", []) => Some(Step::Rotation(q(0), -std::f64::consts::FRAC_PI_4)),
        ("id", []) => gates(&[]),
        _ => None,
    }
}

/// A rotation `exp(-iθ/2 P)`, as `(P, θ)`.
type Rotation = (Pauli, f64);
/// A rotation and the index of the instruction it comes from.
type Located = (Pauli, f64, usize);

/// What [`read`] makes of an instruction: Clifford gates, in order (the
/// first `n` of the array), or one Z rotation. Kept inline rather than in a
/// `Vec`: rotation merging reads every instruction of every pipeline.
enum Step {
    Cliffords([Gate; 6], usize),
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
    let mut rotations = Vec::new();
    for (index, inst) in instructions.iter().enumerate() {
        match read(inst)? {
            Step::Rotation(q, theta) => {
                // Multiples of π/2 are Clifford: S, Z or S†.
                let quarter = theta / FRAC_PI_2;
                if (quarter - quarter.round()).abs() < 1e-12 {
                    let k = quarter.round().rem_euclid(4.0) as usize;
                    apply_all_to_map(&mut m, &[Gate::S(q); 3][..k], &[q]);
                } else {
                    rotations.push((m[2 * q + 1], theta, index));
                }
            }
            Step::Cliffords(gates, n) => apply_all_to_map(&mut m, &gates[..n], &inst.qubits),
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
            Some(theta) => Instruction::gate("rz", inst.qubits, [theta]),
            None => inst,
        })
        .collect()
}

/// How many rotations [`merge_rotations`] looks back past, at most.
const MERGE_LOOKBACK: usize = 256;

/// `C† P C` for the Clifford C whose images `m` holds: the product of the
/// images of P's letters, in the order of the qubits, X before Z on each.
/// Only the qubits in P's support are visited (two at most for the
/// generators [`apply_to_map`] maps), not all of them: on the 32-qubit
/// Quantum Volume circuit, looping over every qubit for every gate was most
/// of rotation merging's time.
fn image(m: &[Pauli], p: &Pauli) -> Pauli {
    let mut out = Pauli {
        phase: p.phase,
        ..Pauli::IDENTITY
    };
    let mut support = p.support();
    while support != 0 {
        let q = support.trailing_zeros() as usize;
        support &= support - 1;
        if (p.x >> q) & 1 == 1 {
            out = out.times(&m[2 * q]);
        }
        if (p.z >> q) & 1 == 1 {
            out = out.times(&m[2 * q + 1]);
        }
    }
    out
}

/// `m ← m ∘ conj(g†)`: appending `g` to the Clifford C makes
/// `(gC)† P (gC) = C† (g† P g) C`. Only the rows of the qubits `g` acts on
/// change.
fn apply_to_map(m: &mut [Pauli], g: Gate) {
    let support: &[usize] = match g {
        Gate::H(q) | Gate::S(q) | Gate::Sdg(q) | Gate::Sx(q) => &[q],
        Gate::Cx(c, t) => &[c, t],
    };
    let inverse = g.inverse();
    // Every new row is computed from the old ones before any is written.
    let mut updates = [(0, Pauli::IDENTITY); 4];
    let mut count = 0;
    for &q in support {
        for (row, generator) in [(2 * q, Pauli::x(q)), (2 * q + 1, Pauli::z(q))] {
            let mut p = generator;
            inverse.conjugate(&mut p);
            updates[count] = (row, image(m, &p));
            count += 1;
        }
    }
    for &(row, p) in &updates[..count] {
        m[row] = p;
    }
}

/// [`apply_to_map`] for each of `gates` in turn, acting on `qubits` (two at
/// most), in one update: `m ← m ∘ conj(g₁†) ∘ … ∘ conj(gₖ†)`. Appending
/// g₁, …, gₖ to C maps each generator P to `C† g₁† ⋯ gₖ† P gₖ ⋯ g₁ C`, so
/// P is conjugated by gₖ† first and g₁† last, and its image read once. The
/// Pauli arithmetic is exact, so the map is the same, bit for bit, as
/// k updates give (`one_update_per_gate_is_the_same_map` checks it). A gate
/// read as several Cliffords (√X as H S H, X as H S S H, CZ as H CX H, Z
/// rotations by π/2 as S, S S or S S S) updates the map once: rotation
/// merging reads every gate of every pipeline, and on the 16-qubit
/// Hamiltonian on the heavy-hex chip with CZ one update per Clifford was
/// most of its time.
fn apply_all_to_map(m: &mut [Pauli], gates: &[Gate], qubits: &[usize]) {
    debug_assert!(qubits.len() <= 2);
    match gates {
        [] => return,
        &[g] => return apply_to_map(m, g),
        _ => {}
    }
    let mut updates = [(0, Pauli::IDENTITY); 4];
    let mut count = 0;
    for &q in qubits {
        for (row, generator) in [(2 * q, Pauli::x(q)), (2 * q + 1, Pauli::z(q))] {
            let mut p = generator;
            for g in gates.iter().rev() {
                g.inverse().conjugate(&mut p);
            }
            updates[count] = (row, image(m, &p));
            count += 1;
        }
    }
    for &(row, p) in &updates[..count] {
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
pub const LOCAL_PENALTY: f64 = 0.01;
/// A penalty that does more than break ties: a gate needing one more local
/// Clifford must remove a tenth of a weight more (a weight in the window's
/// first layer, a third of one in the next) to be chosen. Each weight costs
/// about a CX in the end, a local Clifford about a single-qubit gate, so
/// this gives up a few CXs for fewer single-qubit gates; the pipeline runs
/// trials with both penalties (`network_trials` in `passes.rs`).
pub const FRUGAL_LOCAL_PENALTY: f64 = 0.1;
/// Pairs of qubits tried as places for the next gate, at most.
const MAX_PAIRS: usize = 64;
/// With seed k, at most k times this much random noise is added to each
/// candidate gate's score. Well under a weight, it only reorders gates that
/// do about as well, which is where the greedy synthesis is blind; the seeds
/// spread from gentle to bolder reorderings.
const NOISE: f64 = 0.05;
/// On a device ([`synthesize_on`]): how much the Clifford left counts, in
/// place of [`TABLEAU_WEIGHT`]. There each CX of the network also scatters
/// the frame over the chip, and the rotations still to come are read in it,
/// so a frame kept close to the qubits' own pays for itself. On 12 random
/// Hamiltonians of 6–16 qubits on the heavy-hex chip, with CX: 0.15 left
/// 6298 CXs in all, 0.2 6004, 0.25 6282 (gates 12756, 12032, 11674).
const DEVICE_TABLEAU_WEIGHT: f64 = 0.2;
/// On a device: the gates [`State::best_gate`] may place without a rotation
/// coming out before the first rotation is reduced by force, in place of
/// four per qubit. Wandering costs more on a chip, where moves are long:
/// on the same 12 Hamiltonians, 16 left 6105 CXs, 24 6149, 32 6004, 38
/// (two per qubit of the chip) 6097, and 76 (four per qubit) 6105 with
/// more gates in all.
const DEVICE_PATIENCE: usize = 32;

/// What a Pauli whose qubits are `support` costs on a device: one more than
/// the weight of a minimum spanning tree of `support` under `g`'s distances
/// (its weight, when its qubits are connected on the chip), and at least 1,
/// as the weight is counted at least 1 off a device. Reducing it to one
/// qubit with CXs on coupled pairs takes about that many CXs: a ladder over
/// a Steiner tree of the support would take its edges, which the spanning
/// tree of the support's distances overestimates by at most twice.
fn tree_cost(g: &CouplingGraph, support: u128) -> u32 {
    let mut qubits = [0usize; MAX_QUBITS];
    let mut count = 0;
    let mut rest = support;
    while rest != 0 {
        qubits[count] = rest.trailing_zeros() as usize;
        rest &= rest - 1;
        count += 1;
    }
    if count <= 1 {
        return 1;
    }
    // Prim's algorithm on the complete graph of the support.
    let mut best = [u32::MAX; MAX_QUBITS];
    let mut done = [false; MAX_QUBITS];
    let mut total = 0;
    let mut current = 0;
    done[0] = true;
    for _ in 1..count {
        let mut next = usize::MAX;
        for k in 0..count {
            if done[k] {
                continue;
            }
            best[k] = best[k].min(g.distance(qubits[current], qubits[k]));
            if next == usize::MAX || best[k] < best[next] {
                next = k;
            }
        }
        total += best[next];
        done[next] = true;
        current = next;
    }
    1 + total
}

/// A placement of `instructions`' logical qubits on `graph` for a network
/// synthesized on the chip ([`synthesize_on`]): `layout[q]` is the physical
/// qubit of logical qubit `q`, for `logical` qubits. It minimizes the sum,
/// over the circuit's rotations, of each rotation's spanning-tree cost under
/// the placement ([`tree_cost`], what the synthesis itself scores), by
/// exchanging two physical qubits (an unused one too) while that lowers
/// the sum, from `start` or else from a greedy growth: the logical qubit in
/// most rotations on the most central physical qubit, then, next to the
/// placed qubits, the unplaced one sharing most rotations with a placed
/// neighbour (the placement of Vuorenkoski and Meijer-van de Griend,
/// "Shortest Path in Pauli Forest", arXiv:2605.03545).
///
/// SABRE's layouts are chosen for the circuit as written, CX ladders and
/// all; the network's rotations are what the synthesis on the chip sees.
pub fn placement(
    instructions: &[Instruction],
    logical: usize,
    graph: &CouplingGraph,
    start: Option<&[usize]>,
) -> Option<Vec<usize>> {
    let n = graph.qubits();
    if logical > n || n > MAX_QUBITS {
        return None;
    }
    let (rotations, _) = extract(instructions, logical)?;
    let supports: Vec<u128> = rotations
        .iter()
        .map(|(p, _)| p.support())
        .filter(|s| s.count_ones() > 1)
        .collect();
    if supports.is_empty() {
        return None;
    }
    // at[p]: the logical qubit on physical qubit p, if any.
    let mut at: Vec<Option<usize>> = vec![None; n];
    match start {
        Some(layout) => {
            for (q, &p) in layout.iter().enumerate().take(logical) {
                at[p] = Some(q);
            }
        }
        None => greedy_placement(&supports, logical, graph, &mut at),
    }
    let placed = |at: &[Option<usize>], s: u128| -> u128 {
        let mut out = 0u128;
        for (p, q) in at.iter().enumerate() {
            if let Some(q) = *q
                && s >> q & 1 == 1
            {
                out |= 1 << p;
            }
        }
        out
    };
    let total = |at: &[Option<usize>]| -> u64 {
        supports
            .iter()
            .map(|&s| u64::from(tree_cost(graph, placed(at, s))))
            .sum()
    };
    let mut best = total(&at);
    // First improvement, in a fixed order, until no exchange lowers the sum.
    loop {
        let mut improved = false;
        for a in 0..n {
            for b in a + 1..n {
                if at[a].is_none() && at[b].is_none() {
                    continue;
                }
                at.swap(a, b);
                let cost = total(&at);
                if cost < best {
                    best = cost;
                    improved = true;
                } else {
                    at.swap(a, b);
                }
            }
        }
        if !improved {
            break;
        }
    }
    let mut layout = vec![0; logical];
    for (p, q) in at.iter().enumerate() {
        if let Some(q) = *q {
            layout[q] = p;
        }
    }
    Some(layout)
}

/// The greedy start of [`placement`].
fn greedy_placement(
    supports: &[u128],
    logical: usize,
    graph: &CouplingGraph,
    at: &mut [Option<usize>],
) {
    let n = graph.qubits();
    // shared[a][b]: the rotations on both a and b.
    let mut shared = vec![vec![0u32; logical]; logical];
    let mut uses = vec![0u32; logical];
    for &s in supports {
        for a in 0..logical {
            if s >> a & 1 == 0 {
                continue;
            }
            uses[a] += 1;
            for b in 0..logical {
                if b != a && s >> b & 1 == 1 {
                    shared[a][b] += 1;
                }
            }
        }
    }
    // The most central physical qubit: least total distance to the others,
    // the first on ties.
    let centre = (0..n)
        .min_by_key(|&p| (0..n).map(|r| graph.distance(p, r)).sum::<u32>())
        .expect("a qubit");
    let first = (0..logical)
        .max_by_key(|&q| (uses[q], std::cmp::Reverse(q)))
        .expect("a qubit");
    let mut done = vec![false; logical];
    at[centre] = Some(first);
    done[first] = true;
    for _ in 1..logical {
        // Of the free physical qubits next to placed ones and the unplaced
        // logical qubits, the pair sharing most rotations with that
        // physical qubit's placed neighbours (then most rotations in all).
        let mut choice: Option<((u32, u32), usize, usize)> = None;
        for p in 0..n {
            if at[p].is_some() {
                continue;
            }
            let near: Vec<usize> = graph.neighbours(p).iter().filter_map(|&r| at[r]).collect();
            if near.is_empty() {
                continue;
            }
            for q in (0..logical).filter(|&q| !done[q]) {
                let score = (near.iter().map(|&r| shared[q][r]).sum::<u32>(), uses[q]);
                if choice.as_ref().is_none_or(|(s, _, _)| score > *s) {
                    choice = Some((score, p, q));
                }
            }
        }
        let (_, p, q) = choice.expect("a connected chip");
        at[p] = Some(q);
        done[q] = true;
    }
}

/// Whether the qubits of `set` are connected on `g` through qubits of `set`.
fn connected(g: &CouplingGraph, set: u128) -> bool {
    let start = set.trailing_zeros() as usize;
    let mut seen: u128 = 1 << start;
    let mut stack = vec![start];
    while let Some(u) = stack.pop() {
        for &v in g.neighbours(u) {
            if set >> v & 1 == 1 && seen >> v & 1 == 0 {
                seen |= 1 << v;
                stack.push(v);
            }
        }
    }
    seen == set
}

/// A tree of `g` through qubits of `allowed` only, joining `root` to every
/// qubit of `terminals` (all in `allowed`): each qubit's parent
/// (`usize::MAX` for the root and the qubits outside), and the tree's
/// qubits, each after its parent. Grown from the root by the shortest path
/// to the closest terminal not yet in it, again and again (the heuristic of
/// Takahashi and Matsuyama for Steiner trees, within twice the optimum).
/// None if `allowed` does not connect them.
fn steiner_tree(
    g: &CouplingGraph,
    allowed: u128,
    root: usize,
    terminals: u128,
) -> Option<(Vec<usize>, Vec<usize>)> {
    let n = g.qubits();
    let mut parent = vec![usize::MAX; n];
    let mut previous = vec![usize::MAX; n];
    let mut order = vec![root];
    let mut tree: u128 = 1 << root;
    let mut queue = Vec::with_capacity(n);
    let mut path = Vec::new();
    while terminals & !tree != 0 {
        // Breadth first from the whole tree, to the first terminal met.
        queue.clear();
        queue.extend_from_slice(&order);
        let mut seen = tree;
        let mut found = None;
        let mut head = 0;
        'search: while head < queue.len() {
            let u = queue[head];
            head += 1;
            for &v in g.neighbours(u) {
                if seen >> v & 1 == 1 || allowed >> v & 1 == 0 {
                    continue;
                }
                seen |= 1 << v;
                previous[v] = u;
                queue.push(v);
                if terminals >> v & 1 == 1 {
                    found = Some(v);
                    break 'search;
                }
            }
        }
        let mut v = found?;
        path.clear();
        while tree >> v & 1 == 0 {
            path.push(v);
            v = previous[v];
        }
        for &v in path.iter().rev() {
            parent[v] = previous[v];
            tree |= 1 << v;
            order.push(v);
        }
    }
    Some((parent, order))
}

/// `gate` applied to both Paulis of `ps`, and appended to `gates`.
fn apply(gate: Gate, ps: &mut [Pauli; 2], gates: &mut Vec<Gate>) {
    gate.conjugate(&mut ps[0]);
    gate.conjugate(&mut ps[1]);
    gates.push(gate);
}

/// Gather the letters of `ps[which]` (on the qubits of the tree `parent`,
/// `order` from [`steiner_tree`]) onto the tree's root, children before
/// parents: a letter whose parent has one is cleared by a CX between them
/// once both are Z (Z on the control and the target becomes Z on the
/// target); one whose parent has none is moved there by two (X on the
/// control becomes X on both, and X on both X on the control's partner).
/// So a qubit of the tree outside the Pauli costs two CXs, one in it one.
fn gather(
    ps: &mut [Pauli; 2],
    which: usize,
    parent: &[usize],
    order: &[usize],
    gates: &mut Vec<Gate>,
) {
    for &v in order[1..].iter().rev() {
        let letter = ps[which].letter(v);
        if letter == 0 {
            continue;
        }
        let p = parent[v];
        if ps[which].letter(p) == 0 {
            for gate in to_x(v, letter) {
                apply(gate, ps, gates);
            }
            apply(Gate::Cx(v, p), ps, gates);
            apply(Gate::Cx(p, v), ps, gates);
        } else {
            for gate in to_z(v, letter) {
                apply(gate, ps, gates);
            }
            for gate in to_z(p, ps[which].letter(p)) {
                apply(gate, ps, gates);
            }
            apply(Gate::Cx(v, p), ps, gates);
        }
        debug_assert_eq!(ps[which].letter(v), 0);
    }
}

/// [`to_z`] or [`to_x`].
type LocalTo = fn(usize, u8) -> Vec<Gate>;

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
    synthesize_seeded(instructions, permute, 0, LOCAL_PENALTY)
}

/// [`synthesize`], breaking ties between candidate gates at random when
/// `seed` is not 0, and charging `local_penalty` (in weights) for each
/// single-qubit Clifford a candidate gate needs ([`LOCAL_PENALTY`] in
/// [`synthesize`]). The greedy choice of each CX looks only a few layers
/// ahead, and between gates that do equally well there it can lose a tenth of
/// the CNOTs; a few seeds, run side by side, find the better ones.
pub fn synthesize_seeded(
    instructions: &[Instruction],
    permute: bool,
    seed: u64,
    local_penalty: f64,
) -> Option<(Vec<Instruction>, Vec<usize>)> {
    synthesize_with(instructions, permute, seed, local_penalty, None)
        .map(|mut out| out.swap_remove(0))
}

/// [`synthesize_seeded`] for a device, with relabelling: `instructions` on
/// the physical qubits of `graph`, and every CX that sets up a rotation on a
/// coupled pair. The CXs are chosen by what they do to the rotations'
/// spanning trees on the chip ([`tree_cost`]) rather than their weights, on
/// every coupled pair a rotation of the first layer touches (a CX that
/// grows a rotation towards another qubit of it can lower other rotations'
/// trees), and a rotation whose qubits are apart is reduced by force by
/// moving one of its letters along a shortest path. The Clifford left at
/// the end is synthesized on the chip's pairs too ([`State::clifford_on`]);
/// and, as a second network sharing everything before it, as off a device,
/// to be routed after: of the two, neither needs fewer CXs everywhere.
///
/// Off a device the network is synthesized on logical qubits and then
/// routed, and routing the scattered CXs of a network costs more than the
/// network saves: the 16-qubit benchmark Hamiltonian needs 313 CXs on all
/// pairs, but routed on the heavy-hex chip 1321 (and the direct compilation
/// is kept there). Synthesized on the chip's pairs it needed 1078 with the
/// Clifford routed after (0.4.15), and 1063 with the choice of either.
pub fn synthesize_on(
    instructions: &[Instruction],
    graph: &CouplingGraph,
    seed: u64,
    local_penalty: f64,
) -> Option<Vec<(Vec<Instruction>, Vec<usize>)>> {
    synthesize_with(instructions, true, seed, local_penalty, Some(graph))
}

/// [`synthesize_seeded`], on `graph`'s pairs if given ([`synthesize_on`]).
fn synthesize_with(
    instructions: &[Instruction],
    permute: bool,
    seed: u64,
    local_penalty: f64,
    graph: Option<&CouplingGraph>,
) -> Option<Vec<(Vec<Instruction>, Vec<usize>)>> {
    let used = instructions
        .iter()
        .flat_map(|i| i.qubits.iter().map(|&q| q + 1))
        .max()
        .unwrap_or(0);
    let n = used.max(graph.map_or(0, |g| g.qubits()));
    let (rotations, rows) = extract(instructions, n)?;

    let mut state = State {
        pending: rotations,
        rows,
        out: Vec::new(),
        n,
        rng: (seed != 0).then(|| Rng::new(seed)),
        noise: NOISE * seed as f64,
        local_penalty,
        tables: candidate_tables(),
        scored_buffer: Vec::new(),
        graph,
        trees: Default::default(),
    };
    state.rotations();
    let mut out = Vec::new();
    if let Some(g) = graph.filter(|g| permute && n == g.qubits()) {
        let (rows, network) = (state.rows.clone(), state.out.clone());
        let wires = state.clifford_on(g);
        out.push((std::mem::replace(&mut state.out, network), wires));
        state.rows = rows;
    }
    let wires = state.clifford(permute);
    out.push((state.out, wires));
    Some(out)
}

struct State<'g> {
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
    /// What each single-qubit Clifford of a candidate gate costs.
    local_penalty: f64,
    /// What each candidate gate, σ controlling τ, does to letters
    /// ([`candidate_tables`]), by `[σ − 1][τ − 1]`.
    tables: [[Candidate; 3]; 3],
    /// Reused by [`State::best_gate`] for [`State::scored`]'s Paulis.
    scored_buffer: Vec<Scored>,
    /// The device, if every CX must be on one of its coupled pairs.
    graph: Option<&'g CouplingGraph>,
    /// [`tree_cost`] of each support met so far.
    trees: crate::hashing::FastMap<u128, u32>,
}

/// A Pauli [`State::best_gate`] scores candidates on: how much it counts,
/// its weight and its letters.
struct Scored {
    count: f64,
    weight: u32,
    x: u128,
    z: u128,
}

/// What a candidate controlled-Pauli gate ([`State::controlled_pauli`])
/// does, wherever it is placed: the weight on its two qubits that each of
/// the 16 pairs of letters there (`letter(a) | letter(b) << 2`) has after
/// it, and how many local Cliffords it needs before its CX.
#[derive(Clone, Copy, Debug)]
struct Candidate {
    letters_after: [u32; 16],
    /// The letters on the two qubits after it, as `letter(a) | letter(b) << 2`.
    pairs_after: [u8; 16],
    locals: usize,
}

/// The [`Candidate`] of every σ controlling τ, by `[σ − 1][τ − 1]`. A gate
/// acts on the letters of its own two qubits only, and as on any other two,
/// so these are worked out once on qubits 0 and 1 instead of for every pair
/// and every step: the synthesis scores nine candidates on up to
/// [`MAX_PAIRS`] pairs for each gate it places, and applying their gates to
/// the 16 pairs of letters each time was nearly half of its work.
fn candidate_tables() -> [[Candidate; 3]; 3] {
    let mut tables = [[Candidate {
        letters_after: [0; 16],
        pairs_after: [0; 16],
        locals: 0,
    }; 3]; 3];
    for sigma in 1..=3u8 {
        for tau in 1..=3u8 {
            let gates = State::controlled_pauli(sigma, tau, 0, 1);
            let table = &mut tables[usize::from(sigma - 1)][usize::from(tau - 1)];
            table.locals = gates.len() - 1;
            for (pair, after) in table.letters_after.iter_mut().enumerate() {
                let mut p = Pauli {
                    x: (pair & 1 | (pair >> 2 & 1) << 1) as u128,
                    z: (pair >> 1 & 1 | (pair >> 3 & 1) << 1) as u128,
                    phase: 0,
                };
                gates.iter().for_each(|g| g.conjugate(&mut p));
                *after = p.weight();
                table.pairs_after[pair] = p.letter(0) | (p.letter(1) << 2);
            }
        }
    }
    tables
}

impl State<'_> {
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
                    self.out.push(Instruction::gate("rz", [q], [angle]));
                }
                self.pending.remove(i);
                stuck = 0;
                continue;
            }
            let patience = match self.graph {
                Some(_) => DEVICE_PATIENCE,
                None => 4 * self.n,
            };
            let gates = match stuck < patience {
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

    /// The Paulis whose weights a gate can change: the rotations in the
    /// window, discounted by layer, and the Clifford's images, weighted by
    /// [`TABLEAU_WEIGHT`]. Each as how much it counts, its weight and its
    /// letters, worked out once for every pair [`State::best_gate`] tries.
    fn scored(&self, layers: &[usize], out: &mut Vec<Scored>) {
        let rotations = layers
            .iter()
            .enumerate()
            .map(|(i, &layer)| (DISCOUNT.powi(layer as i32), &self.pending[i].0));
        let weight = match self.graph {
            Some(_) => DEVICE_TABLEAU_WEIGHT,
            None => TABLEAU_WEIGHT,
        };
        let images = self.rows.iter().map(|p| (weight, p));
        out.clear();
        out.extend(rotations.chain(images).map(|(count, p)| Scored {
            count,
            weight: p.weight(),
            x: p.x,
            z: p.z,
        }));
    }

    /// How much each candidate gate on qubits `a` and `b` (σ controlling τ,
    /// at `3(σ − 1) + τ − 1`) changes the cost of the rotations in the
    /// window and of the Clifford left: the weights beyond one of the
    /// [`State::scored`] Paulis that act on `a` or `b`. Only their letters
    /// on `a` and `b` change, and only as those letters say, so each
    /// candidate's table says what becomes of each of the 16 pairs of
    /// letters ([`candidate_tables`]) rather than its gates being applied to
    /// every Pauli. The nine sums are taken in one pass over the Paulis,
    /// each in their order (from −0.0, as `Iterator::sum` does), so each is
    /// the same number as nine passes would give.
    fn changes(tables: &[[Candidate; 3]; 3], scored: &[Scored], a: usize, b: usize) -> [f64; 9] {
        let mask = (1u128 << a) | (1u128 << b);
        let letter = |p: &Scored, q: usize| (((p.x >> q) & 1) | (((p.z >> q) & 1) << 1)) as usize;
        let letters_before = |pair: usize| u32::from(pair & 3 != 0) + u32::from(pair >> 2 != 0);
        let mut sums = [-0.0; 9];
        for p in scored.iter().filter(|p| (p.x | p.z) & mask != 0) {
            let pair = letter(p, a) | (letter(p, b) << 2);
            let unchanged = p.weight - letters_before(pair);
            let before = f64::from(p.weight.max(1));
            for (sum, table) in sums.iter_mut().zip(tables.iter().flatten()) {
                let after = unchanged + table.letters_after[pair];
                *sum += p.count * (f64::from(after.max(1)) - before);
            }
        }
        sums
    }

    /// [`State::changes`] on a device, where a Pauli costs its [`tree_cost`]
    /// instead of its weight. A candidate changes only the letters on `a`
    /// and `b`, so of each Pauli's support only whether `a` and `b` are in
    /// it can change: its four costs are worked out once (and kept in
    /// `trees` by support, as the same supports come back step after step:
    /// on the 16-qubit Hamiltonian this made the synthesis about 6× faster)
    /// and each candidate's table says which applies after it.
    fn changes_on(
        g: &CouplingGraph,
        tables: &[[Candidate; 3]; 3],
        scored: &[Scored],
        a: usize,
        b: usize,
        trees: &mut crate::hashing::FastMap<u128, u32>,
    ) -> [f64; 9] {
        let mask = (1u128 << a) | (1u128 << b);
        let letter = |p: &Scored, q: usize| (((p.x >> q) & 1) | (((p.z >> q) & 1) << 1)) as usize;
        // Which of a and b a pair of letters occupies, as two bits.
        let occupancy =
            |pair: usize| usize::from(pair & 3 != 0) | (usize::from(pair >> 2 != 0) << 1);
        let mut sums = [-0.0; 9];
        for p in scored.iter().filter(|p| (p.x | p.z) & mask != 0) {
            let pair = letter(p, a) | (letter(p, b) << 2);
            let rest = (p.x | p.z) & !mask;
            // The cost for each occupancy of a and b.
            let mut by_occupancy = [0.0; 4];
            for (k, c) in by_occupancy.iter_mut().enumerate() {
                let mut support = rest;
                if k & 1 != 0 {
                    support |= 1 << a;
                }
                if k & 2 != 0 {
                    support |= 1 << b;
                }
                *c = f64::from(
                    *trees
                        .entry(support)
                        .or_insert_with(|| tree_cost(g, support)),
                );
            }
            let before = by_occupancy[occupancy(pair)];
            for (sum, table) in sums.iter_mut().zip(tables.iter().flatten()) {
                let after = by_occupancy[occupancy(usize::from(table.pairs_after[pair]))];
                *sum += p.count * (after - before);
            }
        }
        sums
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
            // On a device, every coupled pair the rotation touches.
            if let Some(g) = self.graph {
                for &a in &qubits {
                    for &b in g.neighbours(a) {
                        let pair = (a.min(b), a.max(b));
                        if !pairs.contains(&pair) {
                            pairs.push(pair);
                        }
                    }
                }
                if pairs.len() >= MAX_PAIRS {
                    break;
                }
                continue;
            }
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
        // (score, change in cost, where) of the best so far, `where` being
        // (σ, τ, a, b). Only the winner's gates are built.
        type Where = (u8, u8, usize, usize);
        let mut best: Option<(f64, f64, Where)> = None;
        let mut scored = std::mem::take(&mut self.scored_buffer);
        self.scored(layers, &mut scored);
        for &(a, b) in &pairs {
            let changes = match self.graph {
                Some(g) => State::changes_on(g, &self.tables, &scored, a, b, &mut self.trees),
                None => State::changes(&self.tables, &scored, a, b),
            };
            for sigma in 1..=3 {
                for tau in 1..=3 {
                    let index = 3 * usize::from(sigma - 1) + usize::from(tau - 1);
                    let table = &self.tables[usize::from(sigma - 1)][usize::from(tau - 1)];
                    let change = changes[index];
                    // Of gates that do as well, those with fewer local
                    // Cliffords before their CX.
                    let mut score = change + self.local_penalty * table.locals as f64;
                    if let Some(rng) = self.rng.as_mut() {
                        score += self.noise * rng.uniform();
                    }
                    if best.as_ref().is_none_or(|(s, _, _)| score < *s - 1e-9) {
                        best = Some((score, change, (sigma, tau, a, b)));
                    }
                }
            }
        }
        self.scored_buffer = scored;
        best.filter(|(_, change, _)| *change < -1e-9)
            .map(|(_, _, (sigma, tau, a, b))| State::controlled_pauli(sigma, tau, a, b))
    }

    /// Lower the weight of rotation `i` by one, whatever it does to the
    /// others; on a device, or bring two of its qubits closer.
    fn reduce_first(&self, i: usize) -> Vec<Gate> {
        let p = self.pending[i].0;
        if let Some(g) = self.graph {
            return State::reduce_on(&p, g);
        }
        let support = p.support();
        let c = support.trailing_zeros() as usize;
        let t = (support & !(1 << c)).trailing_zeros() as usize;
        State::reduction(&p, c, t)
    }

    /// On a device: the closest two qubits of `p`'s support, if coupled,
    /// reduced as [`State::reduction`] does; otherwise the letter of one moved
    /// a step towards the other (Z on `a` and I on `next` becoming I and Z,
    /// two CXs), so that the next call is a step closer.
    fn reduce_on(p: &Pauli, g: &CouplingGraph) -> Vec<Gate> {
        let support: Vec<usize> = (0..g.qubits()).filter(|&q| p.letter(q) != 0).collect();
        let mut closest = (u32::MAX, 0, 0);
        for (k, &a) in support.iter().enumerate() {
            for &b in &support[k + 1..] {
                closest = closest.min((g.distance(a, b), a, b));
            }
        }
        let (d, a, b) = closest;
        if d == 1 {
            return State::reduction(p, a, b);
        }
        let next = g.path(a, b)[1];
        let mut gates = to_z(a, p.letter(a));
        gates.extend([Gate::Cx(next, a), Gate::Cx(a, next)]);
        gates
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
        self.fix_signs(&wires);
        wires
    }

    /// Make each row's pair `+X_w` and `+Z_w` where it is `±X_w` and `±Z_w`
    /// (`w = wires[q]`), with Paulis on `w`.
    fn fix_signs(&mut self, wires: &[usize]) {
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
    }

    /// [`State::clifford`] on a device, with every CX on a coupled pair, after
    /// the architecture-aware synthesis of Clifford tableaus (Winderl et al.,
    /// "Architecture-Aware Synthesis of Stabilizer Circuits from Clifford
    /// Tableaus", arXiv:2309.08972): one row at a time, its pair of images
    /// gathered onto a free qubit `t` along a Steiner tree of the free qubits,
    /// and `t` then taken out of them. `t` is never one whose removal would
    /// split the free qubits, so the rows left stay within one connected
    /// piece of the chip. Of every row and every such `t`, the pair needing
    /// fewest CXs is decoupled first (then fewest gates in all, then the
    /// first). Returns where each row's pair ends.
    ///
    /// The Clifford is what the network leaves of the rotations' frame, on
    /// any qubits; synthesized as off a device and routed, it took about a
    /// tenth of the network's CXs on the chip's Hamiltonians.
    fn clifford_on(&mut self, g: &CouplingGraph) -> Vec<usize> {
        let n = self.n;
        let mut rows: Vec<usize> = (0..n).collect();
        let mut free: u128 = if n == MAX_QUBITS {
            u128::MAX
        } else {
            (1 << n) - 1
        };
        let mut wires = vec![0; n];
        while !rows.is_empty() {
            // (CXs, gates, which row, t, the gates).
            let mut best: Option<(usize, usize, usize, usize, Vec<Gate>)> = None;
            for t in (0..n).filter(|&t| free >> t & 1 == 1) {
                let rest = free & !(1 << t);
                if rest != 0 && !connected(g, rest) {
                    continue;
                }
                for (k, &q) in rows.iter().enumerate() {
                    let pair = [self.rows[2 * q], self.rows[2 * q + 1]];
                    let Some(gates) = State::pair_on(g, pair, t, free) else {
                        continue;
                    };
                    let cxs = gates.iter().filter(|g| matches!(g, Gate::Cx(..))).count();
                    if best
                        .as_ref()
                        .is_none_or(|b| (cxs, gates.len()) < (b.0, b.1))
                    {
                        best = Some((cxs, gates.len(), k, t, gates));
                    }
                }
            }
            let (_, _, k, t, gates) = best.expect("a free qubit that keeps the rest connected");
            for gate in gates {
                self.emit(gate);
            }
            let q = rows.remove(k);
            free &= !(1 << t);
            wires[q] = t;
        }
        self.fix_signs(&wires);
        wires
    }

    /// The gates, CXs on coupled pairs of `free` qubits only, that make
    /// `pair` (a row's images, which anticommute and act on `free` only)
    /// `±X_t` and `±Z_t`, or None if `free` cannot connect them: the cheaper
    /// of [`State::pair_from`] from either image (fewer CXs, then gates).
    fn pair_on(g: &CouplingGraph, pair: [Pauli; 2], t: usize, free: u128) -> Option<Vec<Gate>> {
        let key = |gates: &Vec<Gate>| {
            let cxs = gates.iter().filter(|g| matches!(g, Gate::Cx(..))).count();
            (cxs, gates.len())
        };
        match (
            State::pair_from(g, pair, 0, t, free),
            State::pair_from(g, pair, 1, t, free),
        ) {
            (Some(a), Some(b)) => Some(if key(&b) < key(&a) { b } else { a }),
            (a, b) => a.or(b),
        }
    }

    /// [`State::pair_on`], gathering `pair[first]` first: onto `t` along a
    /// Steiner tree of `free`, where it is made X (Z if `first` is 1); then
    /// the other, which now anticommutes with it on `t`, along a tree of
    /// `free` without `t` onto a neighbour `k` of `t` (the one needing
    /// fewest CXs), from which one CX clears it, keeping `t`'s letter:
    /// CX(k, t) clears Z on `k` and keeps X on `t`, CX(t, k) clears X on
    /// `k` and keeps Z on `t`.
    fn pair_from(
        g: &CouplingGraph,
        pair: [Pauli; 2],
        first: usize,
        t: usize,
        free: u128,
    ) -> Option<Vec<Gate>> {
        let (second, to_first, to_second): (usize, LocalTo, LocalTo) = if first == 0 {
            (1, to_x, to_z)
        } else {
            (0, to_z, to_x)
        };
        let mut ps = pair;
        let mut gates = Vec::new();
        let (parent, order) = steiner_tree(g, free, t, ps[first].support())?;
        gather(&mut ps, first, &parent, &order, &mut gates);
        for gate in to_first(t, ps[first].letter(t)) {
            apply(gate, &mut ps, &mut gates);
        }
        let rest = ps[second].support() & !(1 << t);
        if rest != 0 {
            let allowed = free & !(1 << t);
            let mut best: Option<(usize, [Pauli; 2], Vec<Gate>)> = None;
            for &k in g.neighbours(t).iter().filter(|&&k| allowed >> k & 1 == 1) {
                let Some((parent, order)) = steiner_tree(g, allowed, k, rest) else {
                    continue;
                };
                let (mut ps_k, mut gates_k) = (ps, Vec::new());
                gather(&mut ps_k, second, &parent, &order, &mut gates_k);
                // The letter on k as the second's on t would be cleared by
                // the CX: Z for X on t (CX(k, t)), X for Z on t (CX(t, k)).
                let clearing = if first == 0 { to_z } else { to_x };
                for gate in clearing(k, ps_k[second].letter(k)) {
                    apply(gate, &mut ps_k, &mut gates_k);
                }
                let cx = if first == 0 {
                    Gate::Cx(k, t)
                } else {
                    Gate::Cx(t, k)
                };
                apply(cx, &mut ps_k, &mut gates_k);
                let cxs = gates_k.iter().filter(|g| matches!(g, Gate::Cx(..))).count();
                if best.as_ref().is_none_or(|b| cxs < b.0) {
                    best = Some((cxs, ps_k, gates_k));
                }
            }
            let (_, ps_k, gates_k) = best?;
            ps = ps_k;
            gates.extend(gates_k);
        }
        // The second's letter on t, which anticommutes with the first's, to
        // its own: √X takes Y to Z and keeps X; S† takes Y to X and keeps Z.
        for gate in to_second(t, ps[second].letter(t)) {
            apply(gate, &mut ps, &mut gates);
        }
        debug_assert_eq!((ps[0].x, ps[0].z), (1 << t, 0));
        debug_assert_eq!((ps[1].x, ps[1].z), (0, 1 << t));
        Some(gates)
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
                    "cx" | "cz" | "swap" => Instruction::gate(name, [a, b], []),
                    "rz" => Instruction::gate(name, [a], [rng.range(-3.0, 3.0)]),
                    _ => Instruction::gate(name, [a], []),
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
    fn one_update_per_gate_is_the_same_map() {
        // Every multi-step gate `read` gives, on two qubits of three, in
        // both orders, applied to maps already scrambled by other gates.
        let mut rng = crate::rng::Rng::new(5);
        let pool = [
            Gate::H(0),
            Gate::S(1),
            Gate::Cx(0, 2),
            Gate::Cx(2, 1),
            Gate::H(2),
            Gate::Sdg(0),
        ];
        // The gates `read` gives several Cliffords for, and Z rotations by
        // π/2, π and 3π/2 (S once, twice, three times).
        let mut cases: Vec<(String, Vec<Gate>, Vec<usize>)> = Vec::new();
        for qubits in [[1usize, 2], [2, 0]] {
            for name in ["x", "y", "z", "sx", "sxdg", "cz", "swap"] {
                let arity = if ["cz", "swap"].contains(&name) { 2 } else { 1 };
                let inst = Instruction::gate(name, &qubits[..arity], []);
                let Some(Step::Cliffords(array, n)) = read(&inst) else {
                    panic!("{name} is read as Cliffords");
                };
                assert!(n > 1, "{name}");
                cases.push((name.into(), array[..n].to_vec(), qubits[..arity].to_vec()));
            }
            for k in 2..=3 {
                let q = qubits[0];
                cases.push((format!("S^{k}"), vec![Gate::S(q); k], vec![q]));
            }
        }
        for (name, gates, qubits) in &cases {
            {
                for _ in 0..20 {
                    let mut m: Vec<Pauli> =
                        (0..3).flat_map(|q| [Pauli::x(q), Pauli::z(q)]).collect();
                    for _ in 0..8 {
                        apply_to_map(&mut m, pool[rng.below(pool.len())]);
                    }
                    let mut one = m.clone();
                    apply_all_to_map(&mut one, gates, qubits);
                    for &g in gates {
                        apply_to_map(&mut m, g);
                    }
                    assert_eq!(one, m, "{name} on {qubits:?}");
                }
            }
        }
    }

    #[test]
    fn sx_is_h_s_h() {
        // Every letter on qubit 1 of two, with every phase and a letter on
        // the other qubit that must not change.
        for x in 0..4u128 {
            for z in 0..4u128 {
                for phase in 0..4 {
                    let p = Pauli { x, z, phase };
                    let (mut one, mut three) = (p, p);
                    Gate::Sx(1).conjugate(&mut one);
                    for g in [Gate::H(1), Gate::S(1), Gate::H(1)] {
                        g.conjugate(&mut three);
                    }
                    assert_eq!(one, three, "{p:?}");
                }
            }
        }
    }

    #[test]
    fn candidate_tables_match_the_gates() {
        // The tables `best_gate` scores with, against the gates themselves
        // on two qubits that are not 0 and 1, among others.
        let tables = candidate_tables();
        let (a, b) = (5, 2);
        for sigma in 1..=3u8 {
            for tau in 1..=3u8 {
                let gates = State::controlled_pauli(sigma, tau, a, b);
                let table = &tables[usize::from(sigma - 1)][usize::from(tau - 1)];
                assert_eq!(table.locals, gates.len() - 1);
                for pair in 0..16usize {
                    let mut p = Pauli::IDENTITY;
                    for (q, letter) in [(a, pair & 3), (b, pair >> 2)] {
                        p.x |= ((letter & 1) as u128) << q;
                        p.z |= ((letter >> 1) as u128) << q;
                    }
                    p.x |= 1 << 7;
                    gates.iter().for_each(|g| g.conjugate(&mut p));
                    assert_eq!(table.letters_after[pair] + 1, p.weight());
                }
            }
        }
    }

    #[test]
    fn random_circuits_keep_their_unitary() {
        let mut rng = Rng::new(8);
        for n in [2, 3, 4] {
            for k in 0..60 {
                let circuit = random_circuit(&mut rng, n, 25);
                for permute in [false, true] {
                    // Seeded or not, with either penalty.
                    let seed = (k % 3) as u64;
                    let penalty = [LOCAL_PENALTY, FRUGAL_LOCAL_PENALTY][k % 2];
                    let (mut out, wires) = synthesize_seeded(&circuit, permute, seed, penalty)
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
                            out.push(Instruction::gate("swap", [w, q], []));
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
    fn tree_costs_are_spanning_trees_on_the_chip() {
        // 0 – 1 – 2 – 3 – 4 – 5, a line.
        let line = CouplingGraph::line(6);
        let support = |qubits: &[usize]| qubits.iter().fold(0u128, |s, &q| s | 1 << q);
        assert_eq!(tree_cost(&line, 0), 1);
        assert_eq!(tree_cost(&line, support(&[3])), 1);
        // Neighbours: the weight.
        assert_eq!(tree_cost(&line, support(&[1, 2, 3])), 3);
        // Apart: one more for each qubit between.
        assert_eq!(tree_cost(&line, support(&[0, 5])), 6);
        assert_eq!(tree_cost(&line, support(&[0, 2, 5])), 6);
        // A 2 × 3 grid, 0 1 2 over 3 4 5: 0, 2 and 4 are 2 apart in pairs,
        // so the spanning tree weighs 4 (a Steiner tree through 1 would
        // weigh 3: the cost overestimates, by at most twice).
        let grid = CouplingGraph::grid(2, 3);
        assert_eq!(tree_cost(&grid, support(&[0, 2, 4])), 5);
        assert_eq!(tree_cost(&grid, support(&[0, 1, 2, 4])), 4);
    }

    /// `out` followed by SWAPs moving each state from `wires[q]` back to `q`.
    fn unpermuted(mut out: Vec<Instruction>, wires: &[usize]) -> Vec<Instruction> {
        let mut at = wires.to_vec();
        for q in 0..at.len() {
            let w = at[q];
            if w != q {
                out.push(Instruction::gate("swap", [w, q], []));
                let p = at.iter().position(|&x| x == q).expect("a permutation");
                at[p] = w;
                at[q] = q;
            }
        }
        out
    }

    /// A line, a 2 × 3 grid, and a tree with cut vertices everywhere but its
    /// leaves (0 – 1 – 3 – 4 – 5, with 2 on 1 and 6 on 4), where the Clifford
    /// on the chip must end its rows on leaves first.
    fn small_chips() -> [CouplingGraph; 3] {
        let tree = [(0, 1), (1, 2), (1, 3), (3, 4), (4, 5), (4, 6)];
        [
            CouplingGraph::line(5),
            CouplingGraph::grid(2, 3),
            CouplingGraph::new(7, &tree).expect("a tree"),
        ]
    }

    #[test]
    fn networks_on_a_chip_keep_their_unitary_and_set_up_rotations_on_its_pairs() {
        let mut rng = Rng::new(12);
        for graph in &small_chips() {
            let n = graph.qubits();
            for k in 0..12 {
                // Pauli exponentials on the first qubits, idle ones beyond;
                // or Clifford gates on any pairs, with a few rotations.
                let circuit = if k % 4 == 3 {
                    random_circuit(&mut rng, n, 40)
                } else {
                    crate::testing::pauli_exponentials(&mut rng, n - k % 2, 6)
                };
                let seed = (k % 3) as u64;
                let penalty = [LOCAL_PENALTY, FRUGAL_LOCAL_PENALTY][k % 2];
                let networks = synthesize_on(&circuit, graph, seed, penalty).expect("a network");
                assert_eq!(networks.len(), 2);
                for (kind, (out, wires)) in networks.into_iter().enumerate() {
                    // With the Clifford synthesized on the chip, every CX is on
                    // a coupled pair; with it synthesized as off a device,
                    // every CX up to the last rotation.
                    let last = match kind {
                        0 => out.len(),
                        _ => out.iter().rposition(|i| i.name == "rz").unwrap_or(0),
                    };
                    for inst in out[..last].iter().filter(|i| i.name == "cx") {
                        assert!(graph.coupled(inst.qubits[0], inst.qubits[1]), "{inst:?}");
                    }
                    assert!(dense_equal_up_to_phase(
                        &circuit_unitary(&unpermuted(out, &wires), n),
                        &circuit_unitary(&circuit, n),
                        1e-8
                    ));
                }
            }
        }
    }

    /// The sum [`placement`] lowers: each rotation's tree cost, placed.
    fn placed_cost(
        circuit: &[Instruction],
        logical: usize,
        g: &CouplingGraph,
        layout: &[usize],
    ) -> u32 {
        let (rotations, _) = extract(circuit, logical).expect("rotations");
        rotations
            .iter()
            .map(|(p, _)| p.support())
            .filter(|s| s.count_ones() > 1)
            .map(|s| {
                let placed = (0..logical)
                    .filter(|&q| s >> q & 1 == 1)
                    .fold(0u128, |out, q| out | 1 << layout[q]);
                tree_cost(g, placed)
            })
            .sum()
    }

    #[test]
    fn placements_are_layouts_no_costlier_than_their_start() {
        let mut rng = Rng::new(31);
        let mut better = 0;
        for graph in &small_chips() {
            let n = graph.qubits();
            for k in 0..8 {
                let logical = n - k % 2;
                let circuit = crate::testing::pauli_exponentials(&mut rng, logical, 6);
                let start: Vec<usize> = (0..logical).rev().collect();
                for from in [None, Some(&start[..])] {
                    let layout = placement(&circuit, logical, graph, from).expect("a placement");
                    assert_eq!(layout.len(), logical);
                    let mut seen = vec![false; n];
                    for &p in &layout {
                        assert!(p < n && !seen[p], "{layout:?}");
                        seen[p] = true;
                    }
                    if from.is_some() {
                        let (after, before) = (
                            placed_cost(&circuit, logical, graph, &layout),
                            placed_cost(&circuit, logical, graph, &start),
                        );
                        assert!(after <= before, "{after} > {before}");
                        better += usize::from(after < before);
                    }
                }
            }
        }
        // Not vacuous: the exchanges improve most starts.
        assert!(better >= 12, "{better}");
    }

    #[test]
    fn placements_put_a_chain_of_rotations_on_a_line() {
        // ZZ rotations on the neighbours of a shuffled chain of 6 qubits:
        // on a line, the best placement follows the chain, and every
        // rotation's tree costs 2.
        let chain = [3, 0, 5, 1, 4, 2];
        let mut circuit = Vec::new();
        for w in chain.windows(2) {
            circuit.push(Instruction::gate("cx", [w[0], w[1]], []));
            circuit.push(Instruction::gate("rz", [w[1]], [0.3]));
            circuit.push(Instruction::gate("cx", [w[0], w[1]], []));
        }
        let line = CouplingGraph::line(6);
        for from in [None, Some(&[0, 1, 2, 3, 4, 5][..])] {
            let layout = placement(&circuit, 6, &line, from).expect("a placement");
            assert_eq!(placed_cost(&circuit, 6, &line, &layout), 10, "{layout:?}");
            for w in chain.windows(2) {
                assert!(line.coupled(layout[w[0]], layout[w[1]]), "{layout:?}");
            }
        }
    }

    #[test]
    fn cliffords_on_a_chip_need_coupled_pairs_only() {
        // Random Clifford circuits alone (no rotation to set up), so the
        // network is the Clifford synthesis only; on the tree its rows must
        // end on leaves, never on a qubit that would split the free ones.
        let mut rng = Rng::new(31);
        let mut fewer = 0;
        for graph in &small_chips() {
            let n = graph.qubits();
            for _ in 0..20 {
                let circuit: Vec<Instruction> = random_circuit(&mut rng, n, 60)
                    .into_iter()
                    .filter(|i| !["rz", "t"].contains(&i.name.as_str()))
                    .collect();
                let networks = synthesize_on(&circuit, graph, 0, LOCAL_PENALTY).expect("a network");
                let cxs = |out: &[Instruction]| out.iter().filter(|i| i.name == "cx").count();
                // The other, routed from the identity layout: its CXs plus
                // the CXs of its SWAPs.
                let identity: Vec<usize> = (0..n).collect();
                let routed = crate::routing::route(
                    &networks[1].0,
                    n,
                    graph,
                    &identity,
                    &crate::routing::SabreOptions::default(),
                )
                .expect("a routing");
                let routed_cxs = cxs(&routed.instructions) + routed.swap_cost;
                let (out, wires) = &networks[0];
                assert!(
                    out.iter()
                        .filter(|i| i.name == "cx")
                        .all(|i| graph.coupled(i.qubits[0], i.qubits[1]))
                );
                if cxs(out) < routed_cxs {
                    fewer += 1;
                }
                assert!(dense_equal_up_to_phase(
                    &circuit_unitary(&unpermuted(out.clone(), wires), n),
                    &circuit_unitary(&circuit, n),
                    1e-8
                ));
            }
        }
        // Not vacuous, and the point of it: on these chips the Clifford
        // synthesized on them mostly needs fewer CXs than the other routed.
        assert!(fewer > 30, "{fewer} of 60");
    }

    #[test]
    fn steiner_trees_join_their_terminals_through_allowed_qubits() {
        let grid = CouplingGraph::grid(3, 3);
        let support = |qs: &[usize]| qs.iter().fold(0u128, |s, &q| s | 1 << q);
        // Corners 0, 2 and 8 of a 3 × 3 grid from the centre, avoiding 1.
        let allowed = support(&[0, 2, 3, 4, 5, 6, 7, 8]);
        let (parent, order) = steiner_tree(&grid, allowed, 4, support(&[0, 2, 8])).expect("a tree");
        assert_eq!(order[0], 4);
        for &v in &order[1..] {
            assert!(grid.coupled(v, parent[v]));
            assert!(
                order.iter().position(|&p| p == parent[v]) < order.iter().position(|&p| p == v)
            );
        }
        assert!(order.iter().all(|&v| allowed >> v & 1 == 1));
        assert!([0, 2, 8].iter().all(|t| order.contains(t)));
        // 4–3–0, 4–5–2, 5–8: six qubits.
        assert_eq!(order.len(), 6);
        // Cut off: 2 cannot be reached without 1 and 5.
        assert!(steiner_tree(&grid, support(&[0, 2, 3, 4]), 4, support(&[2])).is_none());
        assert!(connected(&grid, support(&[0, 1, 2])));
        assert!(!connected(&grid, support(&[0, 2, 4])));
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
        let g = |name: &str, q: Vec<usize>| Instruction::gate(name, q, []);
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
