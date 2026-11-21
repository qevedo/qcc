//! Pairs of single-qubit gates around two-qubit gates, chosen wire by wire:
//! Paulis, and rotations about the axis a gate commutes with.
//!
//! If a single-qubit gate `D` on one qubit commutes with a two-qubit gate
//! `G` (X with a CX's target, Z with its control or with either qubit of a
//! CZ, and so any rotation about those axes), then `G = D† G D`: a `D` can
//! be put on that qubit just before `G` and a `D†` just after it without
//! changing the circuit. Each merges into the single-qubit run on its side,
//! so the run before `G` becomes `D R` and the one after `R' D†`. Doing so
//! around several gates along a wire carries a gate from one run to
//! another, however many such gates lie between: an X left behind by a
//! block's synthesis, alone between two CXs it is the target of, can travel
//! to a run further on that absorbs it for free (in {rz, sx, x},
//! `X·Rz(b)·SX·Rz(a)` is `Rz(π − b)·SX·Rz(π + a)`, still three gates), or
//! meet another X and cancel; the last Z rotation of a run can travel along
//! a qubit through the controls of several CXs, or through CZs, into a run
//! that starts with one, and the two become one gate.
//!
//! Each pair leaves the circuit's unitary as it was, and pairs on different
//! wires do not interact, so every wire is solved on its own and exactly
//! over the choices it is offered: the wire is a chain of runs separated by
//! the gates on it, the choice at each gate is the `D` of the pair there,
//! and the cost of a run depends only on the choices at the gates on its two
//! sides. A dynamic programme over the chain finds the choices that make the
//! runs cheapest in all; a wire changes only if that is strictly cheaper
//! than leaving it alone, and only its changed runs are resynthesized.
//!
//! The choices at a gate: the identity always; where a whole axis commutes
//! with the gate on that qubit, the Pauli of that axis and the rotations
//! about it that would remove an outer rotation (in the same axis) of one of
//! the nearest runs on either side that the rotation can reach through
//! gates of the same axis ([`WINDOW`] of them each way); elsewhere the
//! Paulis that commute with it (an ECR's). The Paulis alone found what
//! 0.4.13 shipped; the rotations are a superset of them, so a wire never
//! comes out costlier than it would with Paulis alone.
//!
//! Block resynthesis and the shifts through blocks (`symmetric_blocks.rs`)
//! look at one block and the runs just around it; a gate that pays only two
//! or more gates further along the wire is out of their sight. On QCC's
//! final circuits these are mostly the X gates the syntheses leave (the
//! 32-qubit Heisenberg chain had 66), and Z rotations left on either side
//! of a chain of CZs or CX controls.

use crate::gates::{h, matrix_1q, matrix_2q, pauli, rotation};
use crate::hashing::FastMap;
use crate::ir::{Instruction, Name};
use crate::linalg::{M2, M4, equal_up_to_phase4, kron};
use crate::one_qubit::{Basis1q, TOL, count_1q, euler_zyz, is_zero, synthesize_1q, wrap};
use std::f64::consts::{FRAC_PI_2, PI};

/// How many runs on each side of a gate offer it their outer rotations:
/// the nearest non-empty ones its axis reaches. A rotation offered by a run
/// further away is still found if every gate between is offered it, which
/// the nearer gates are.
const WINDOW: usize = 2;

/// The identity and the three Paulis, indexed 0 (I), 1 (X), 2 (Y), 3 (Z).
fn paulis() -> [M2; 4] {
    [M2::IDENTITY, pauli(0), pauli(1), pauli(2)]
}

/// A run of single-qubit gates on one qubit, or any other instruction.
enum Item {
    Run {
        qubit: usize,
        m: M2,
        gates: Vec<Instruction>,
    },
    Op(Instruction),
}

/// A place on a wire where a run is or could be: between two operations on
/// the qubit, or before the first, or after the last.
struct Slot {
    /// The run there, by item index.
    run: Option<usize>,
    /// The operation that ends the slot, by item index; None at the end.
    until: Option<usize>,
}

/// What may be put around an operation on one of its qubits.
#[derive(Clone, Copy, PartialEq, Debug)]
enum Side {
    /// Every rotation about this axis (0 for X, 2 for Z, as in
    /// [`rotation`]) commutes with the gate there.
    Axis(usize),
    /// Only these Paulis do (bit `k` for Pauli `k` of [`paulis`]); `1`,
    /// the identity alone, for anything that is not a two-qubit gate.
    Paulis(u8),
}

/// The sides of [`sides`] by gate name and parameters: a circuit has a few
/// kinds of two-qubit gate (a CX for most), each tried once.
type Known = FastMap<(Name, [u64; 3]), [Side; 2]>;

/// For each qubit of `inst`, what commutes with it there, up to a global
/// phase. Only the identity for anything but a two-qubit gate with a known
/// matrix: a measurement, a reset, a barrier or a gate on three qubits
/// stops everything.
fn sides(inst: &Instruction, paulis: &[M2; 4], known: &mut Known) -> Vec<Side> {
    if inst.qubits.len() != 2 || !inst.is_gate() {
        return vec![Side::Paulis(1); inst.qubits.len()];
    }
    let mut bits = [0u64; 3];
    for (b, p) in bits.iter_mut().zip(inst.params.as_slice()) {
        *b = p.to_bits();
    }
    let sides = known.entry((inst.name.clone(), bits)).or_insert_with(|| {
        let mut sides = [Side::Paulis(1); 2];
        let Some(g) = matrix_2q(&inst.name, &inst.params) else {
            return sides;
        };
        let commutes = |k: usize, d: &M2| {
            let on: M4 = if k == 0 {
                kron(d, &M2::IDENTITY)
            } else {
                kron(&M2::IDENTITY, d)
            };
            equal_up_to_phase4(&(on.adjoint() * g * on), &g, 1e-9)
        };
        for (k, side) in sides.iter_mut().enumerate() {
            // A rotation by a generic angle commutes only if its axis does.
            if let Some(axis) = [2, 0].into_iter().find(|&a| commutes(k, &rotation(a, 0.7))) {
                *side = Side::Axis(axis);
                continue;
            }
            let mut mask = 1u8;
            for (p, pm) in paulis.iter().enumerate().skip(1) {
                if commutes(k, pm) {
                    mask |= 1 << p;
                }
            }
            *side = Side::Paulis(mask);
        }
        sides
    });
    sides.to_vec()
}

/// `m` seen in the Z basis when `axis` is X (`H m H`), as itself for Z.
fn in_z(m: &M2, axis: usize) -> M2 {
    if axis == 2 { *m } else { h() * *m * h() }
}

/// Angles `a` for which `m · R(a)†` (`first`: the run after the gate,
/// where `R(a)†` comes first in time) or `R(a) · m` (the run before it)
/// loses an outer rotation about the axis, given the ZYZ angles of `m` in
/// that axis, `Rz(φ) Ry(θ) Rz(λ)`: `λ` and `−φ`, each also shifted by π
/// (the same unitary with `−θ`) and by −π/2 (the outer angles of the
/// {rz, sx} form at θ = π/2, `Rz(φ + π/2) SX Rz(λ − π/2)`); a run that is a
/// rotation about the axis alone offers all of it.
fn edge_angles((theta, phi, lam): (f64, f64, f64), first: bool, out: &mut Vec<f64>) {
    if is_zero(theta) {
        out.push(if first { phi + lam } else { -(phi + lam) });
        return;
    }
    let a = if first { lam } else { -phi };
    out.extend([a, a + PI, a - FRAC_PI_2]);
}

/// The choices at every gate of one wire, one after the other (see the
/// module's description): those of gate `g` run from `start[g]` to
/// `start[g + 1]` in `options`. The identity is always first, then the Paulis, then, with
/// rotations, rotations by angles neither zero nor π. Each comes with its
/// angle about the gate's axis where the gate has one (zero for the
/// identity, π for the Pauli), NaN elsewhere.
#[derive(Default)]
struct Choices {
    options: Vec<(M2, f64)>,
    start: Vec<usize>,
}

impl Choices {
    fn at(&self, g: usize) -> &[(M2, f64)] {
        &self.options[self.start[g]..self.start[g + 1]]
    }

    /// Fill in the choices of a wire: `sides[g]` is what commutes with gate
    /// `g` there, slot `s` lies between gate `s - 1` and gate `s`, and
    /// `euler(s, axis)` gives the ZYZ angles of slot `s`'s run in that axis,
    /// or None for an empty slot.
    fn fill(
        &mut self,
        rotations: bool,
        sides: &[Side],
        slots: usize,
        euler: &mut dyn FnMut(usize, usize) -> Option<(f64, f64, f64)>,
        ps: &[M2; 4],
    ) {
        self.options.clear();
        self.start.clear();
        let mut angles = Vec::new();
        let mut kept: Vec<f64> = Vec::new();
        for (g, side) in sides.iter().enumerate() {
            self.start.push(self.options.len());
            let axis = match *side {
                Side::Paulis(mask) => {
                    let allowed = (0..4).filter(|p| mask & (1 << p) != 0);
                    self.options.extend(allowed.map(|p| (ps[p], f64::NAN)));
                    continue;
                }
                Side::Axis(axis) => axis,
            };
            self.options
                .extend([(M2::IDENTITY, 0.0), (pauli(axis), PI)]);
            if !rotations {
                continue;
            }
            angles.clear();
            // Forward: slot g + 1 is right after gate g; slot s is reached
            // while the gates before it, from g on, share the axis.
            let mut found = 0;
            for s in g + 1..slots {
                if let Some(e) = euler(s, axis) {
                    edge_angles(e, true, &mut angles);
                    found += 1;
                }
                if found == WINDOW || s >= sides.len() || sides[s] != Side::Axis(axis) {
                    break;
                }
            }
            // Back: slot g is right before gate g.
            let mut found = 0;
            for s in (0..=g).rev() {
                if let Some(e) = euler(s, axis) {
                    edge_angles(e, false, &mut angles);
                    found += 1;
                }
                if found == WINDOW || s == 0 || sides[s - 1] != Side::Axis(axis) {
                    break;
                }
            }
            kept.clear();
            for &a in &angles {
                let a = wrap(a);
                // Both in (−π, π], away from ±π: equal modulo 2π only if
                // close as numbers.
                if is_zero(a) || is_zero(a - PI) || kept.iter().any(|&k| (k - a).abs() < TOL) {
                    continue;
                }
                kept.push(a);
                self.options.push((rotation(axis, a), a));
            }
        }
        self.start.push(self.options.len());
    }
}

/// Insert pairs of Paulis around the two-qubit gates where that makes the
/// single-qubit runs cheaper in `bases` (see the module's description).
/// Changed runs are resynthesized, new ones created where a wire had none;
/// every other run is left as it was.
pub fn insert_pauli_pairs(instructions: Vec<Instruction>, bases: &[Basis1q]) -> Vec<Instruction> {
    insert_pairs(instructions, bases, false, true)
}

/// As [`insert_pauli_pairs`], with rotations about the axis that commutes
/// with each gate offered as well as its Pauli. Never costlier than
/// [`insert_pauli_pairs`] on the same circuit, nor than the circuit itself.
pub fn insert_rotation_pairs(
    instructions: Vec<Instruction>,
    bases: &[Basis1q],
) -> Vec<Instruction> {
    insert_pairs(instructions, bases, true, true)
}

/// The pairs, with rotations or Paulis only; `split` allows pricing a run's
/// two sides apart where its count allows (off only to test that it
/// changes nothing).
fn insert_pairs(
    instructions: Vec<Instruction>,
    bases: &[Basis1q],
    rotations: bool,
    split: bool,
) -> Vec<Instruction> {
    if bases.is_empty() {
        return instructions;
    }
    let ps = paulis();
    let mut known = Known::default();
    let mut items: Vec<Item> = Vec::with_capacity(instructions.len());
    // Per qubit: its slots so far, what commutes with the operation that
    // ends each closed slot, and whether its last slot holds a run that is
    // still open.
    let mut slots: Vec<Vec<Slot>> = Vec::new();
    let mut allowed: Vec<Vec<Side>> = Vec::new();
    let grow = |slots: &mut Vec<Vec<Slot>>, allowed: &mut Vec<Vec<Side>>, q: usize| {
        while slots.len() <= q {
            slots.push(vec![Slot {
                run: None,
                until: None,
            }]);
            allowed.push(Vec::new());
        }
    };
    for inst in instructions {
        let matrix = match inst.qubits.as_slice() {
            [_] if inst.is_gate() => matrix_1q(&inst.name, &inst.params),
            _ => None,
        };
        if let Some(g) = matrix {
            let q = inst.qubits[0];
            grow(&mut slots, &mut allowed, q);
            let slot = slots[q].last_mut().expect("a slot per qubit");
            match slot.run {
                Some(r) => {
                    if let Item::Run { m, gates, .. } = &mut items[r] {
                        *m = g * *m;
                        gates.push(inst);
                    }
                }
                None => {
                    slot.run = Some(items.len());
                    items.push(Item::Run {
                        qubit: q,
                        m: g,
                        gates: vec![inst],
                    });
                }
            }
            continue;
        }
        let at = items.len();
        for (&q, side) in inst.qubits.iter().zip(sides(&inst, &ps, &mut known)) {
            grow(&mut slots, &mut allowed, q);
            slots[q].last_mut().expect("a slot per qubit").until = Some(at);
            slots[q].push(Slot {
                run: None,
                until: None,
            });
            allowed[q].push(side);
        }
        items.push(Item::Op(inst));
    }

    // New matrices for runs that change, by item index, and runs to create,
    // as (item index to put them before, or None for the end; qubit; matrix).
    let mut changed: Vec<Option<M2>> = vec![None; items.len()];
    let mut created: Vec<(Option<usize>, usize, M2)> = Vec::new();
    // Buffers reused from wire to wire.
    let mut choices = Choices::default();
    let mut eulers: Vec<[Option<(f64, f64, f64)>; 2]> = Vec::new();
    let (mut best, mut next, mut from) = (Vec::new(), Vec::new(), Vec::new());
    let xyx = bases.iter().any(|b| matches!(b, Basis1q::Xyx));
    for (q, (wire, sides)) in slots.iter().zip(&allowed).enumerate() {
        if sides.iter().all(|&s| s == Side::Paulis(1)) {
            continue;
        }
        let run_of = |slot: &Slot| match slot.run.map(|r| &items[r]) {
            Some(Item::Run { m, .. }) => Some(*m),
            _ => None,
        };
        let len_of = |slot: &Slot| match slot.run.map(|r| &items[r]) {
            Some(Item::Run { gates, .. }) => gates.len(),
            _ => 0,
        };
        // The ZYZ angles of each run in each axis (Z, and X through `H m
        // H`), worked out when first asked for: a run offers its outer
        // rotations to several gates around it.
        eulers.clear();
        eulers.resize(wire.len(), [None; 2]);
        let mut euler = |s: usize, axis: usize| -> Option<(f64, f64, f64)> {
            let m = run_of(&wire[s])?;
            let k = usize::from(axis == 0);
            Some(*eulers[s][k].get_or_insert_with(|| euler_zyz(&in_z(&m, axis))))
        };
        choices.fill(rotations, sides, wire.len(), &mut euler, &ps);
        // Where the gates on both sides of a slot (or the wire's end) leave
        // only Z rotations to choose from, `Rz(α) · Rz(φ) Ry(θ) Rz(λ) ·
        // Rz(−β)` is `Rz(φ + α) Ry(θ) Rz(λ − β)`: the run's ZYZ angles,
        // worked out once, give every pair's count without a product or
        // another decomposition. Most slots of a CZ circuit, and the
        // controls of a CX one, are of this kind.
        let z_side = |g: Option<usize>| g.is_none_or(|g| sides[g] == Side::Axis(2));
        let by_angles: Vec<Option<(f64, f64, f64)>> = (0..wire.len())
            .map(|s| {
                let both = rotations
                    && z_side(s.checked_sub(1))
                    && z_side((s < sides.len()).then_some(s))
                    && !xyx;
                both.then(|| euler(s, 2).unwrap_or((0.0, 0.0, 0.0)))
            })
            .collect();
        let choices = &choices;
        // The choice `c` at gate `g`, or the identity beyond either end.
        let choice = |g: Option<usize>, c: usize| match g {
            Some(g) if g < sides.len() => choices.at(g)[c],
            _ => (M2::IDENTITY, 0.0),
        };
        // The cost of slot `s` with `before` (by index into the choices of
        // gate s - 1) and `after` (into gate s's): as it is when both are
        // the identity, else as resynthesized.
        let slot_cost = |s: usize, before: usize, after: usize| -> usize {
            if before == 0 && after == 0 {
                return len_of(&wire[s]);
            }
            let (b, beta) = choice(s.checked_sub(1), before);
            let (a, alpha) = choice(Some(s), after);
            if let Some((theta, phi, lam)) = by_angles[s] {
                let angles = (theta, phi + alpha, lam - beta);
                // Not read by the bases this path admits (all but Xyx).
                let unused = M2::IDENTITY;
                return bases
                    .iter()
                    .map(|basis| basis.count_from_angles(&unused, angles))
                    .min()
                    .unwrap_or(0);
            }
            match run_of(&wire[s]) {
                Some(m) => count_1q(&(a * m * b.adjoint()), bases),
                None => count_1q(&(a * b.adjoint()), bases),
            }
        };
        let unchanged: usize = (0..wire.len()).map(|s| slot_cost(s, 0, 0)).sum();
        // best[c]: the cheapest cost of the slots up to the current gate,
        // with choice c there; from[start[g] + c]: the choice at the gate
        // before gate g.
        best.clear();
        best.push(0);
        from.clear();
        from.resize(choices.options.len(), 0);
        for s in 0..sides.len() {
            let here = choices.start[s];
            let count = choices.start[s + 1] - here;
            next.clear();
            next.resize(count, usize::MAX);
            // A run between two Z-axis gates whose count is a sum of a term
            // for each side (`Basis1q::split_count`): the cheapest choice
            // before is the same for every choice here, so it is found once
            // (K prices instead of K²), with the same ties as below. The
            // identity on both sides still costs the run as it stands.
            let sides_apart = match (by_angles[s], bases) {
                (Some(angles), [basis]) if split => {
                    basis.split_count(angles.0).map(|d| (angles, d))
                }
                _ => None,
            };
            if let Some(((_, phi, lam), (base, d_lam, d_phi))) = sides_apart {
                // (total, choice) over every choice before, and over all but
                // the identity.
                let (mut all, mut rest) = ((usize::MAX, 0), (usize::MAX, 0));
                for (before, &so_far) in best.iter().enumerate() {
                    let beta = choice(s.checked_sub(1), before).1;
                    let total = so_far + usize::from(!is_zero(lam - beta + d_lam));
                    if total < all.0 {
                        all = (total, before);
                    }
                    if before > 0 && total < rest.0 {
                        rest = (total, before);
                    }
                }
                for (after, slot) in next.iter_mut().enumerate() {
                    let alpha = choices.at(s)[after].1;
                    let here_cost = base + usize::from(!is_zero(phi + alpha + d_phi));
                    let (total, back) = if after > 0 {
                        (all.0 + here_cost, all.1)
                    } else {
                        let unchanged = best[0] + len_of(&wire[s]);
                        match rest.0.checked_add(here_cost) {
                            Some(t) if t < unchanged => (t, rest.1),
                            _ => (unchanged, 0),
                        }
                    };
                    *slot = total;
                    from[here + after] = back;
                }
                std::mem::swap(&mut best, &mut next);
                continue;
            }
            for (after, slot) in next.iter_mut().enumerate() {
                // The identity first and strictly cheaper only: ties keep
                // the earlier choice. A slot costs nothing or more, so a
                // choice before that has already cost as much as the best
                // total found cannot win, and its slot is not priced.
                for (before, &so_far) in best.iter().enumerate() {
                    if so_far >= *slot {
                        continue;
                    }
                    let total = so_far + slot_cost(s, before, after);
                    if total < *slot {
                        *slot = total;
                        from[here + after] = before;
                    }
                }
            }
            std::mem::swap(&mut best, &mut next);
        }
        let last = wire.len() - 1;
        let (total, mut before) = (0..best.len())
            .map(|c| (best[c] + slot_cost(last, c, 0), c))
            .min()
            .expect("the identity is always a choice");
        if total >= unchanged {
            continue;
        }
        // Walk back from the end: slot s lies between gate s - 1 and gate s.
        let mut after = 0;
        for s in (0..wire.len()).rev() {
            if before != 0 || after != 0 {
                let (b, _) = choice(s.checked_sub(1), before);
                let (a, _) = choice(Some(s), after);
                let m = a * run_of(&wire[s]).unwrap_or(M2::IDENTITY) * b.adjoint();
                match wire[s].run {
                    Some(r) => changed[r] = Some(m),
                    None => created.push((wire[s].until, q, m)),
                }
            }
            if s > 0 {
                after = before;
                before = from[choices.start[s - 1] + before];
            }
        }
    }
    if created.is_empty() && changed.iter().all(Option::is_none) {
        return items
            .into_iter()
            .flat_map(|item| match item {
                Item::Run { gates, .. } => gates,
                Item::Op(inst) => vec![inst],
            })
            .collect();
    }

    // Runs created before the same operation, and at the end, go in the
    // order of their qubits, so that the output does not depend on more
    // than the input.
    created.sort_by_key(|&(at, q, _)| (at.is_none(), at, q));
    let mut created = created.into_iter().peekable();
    let emit = |qubit: usize, m: &M2| -> Vec<Instruction> {
        synthesize_1q(m, bases)
            .into_iter()
            .map(|op| Instruction::gate(&op.name, [qubit], op.params))
            .collect()
    };
    let mut out = Vec::new();
    for (i, item) in items.into_iter().enumerate() {
        while let Some((_, q, m)) = created.next_if(|&(at, _, _)| at == Some(i)) {
            out.extend(emit(q, &m));
        }
        match item {
            Item::Run { qubit, gates, .. } => match changed[i] {
                Some(m) => out.extend(emit(qubit, &m)),
                None => out.extend(gates),
            },
            Item::Op(inst) => out.push(inst),
        }
    }
    for (_, q, m) in created {
        out.extend(emit(q, &m));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ir::cost;
    use crate::lower::Lowering;
    use crate::one_qubit::bases_for;
    use crate::testing::{Rng, circuit_unitary, dense_equal_up_to_phase};
    use std::collections::HashSet;
    use std::f64::consts::FRAC_PI_2;

    fn bases(names: &[&str]) -> Vec<Basis1q> {
        bases_for(&names.iter().map(|s| s.to_string()).collect::<HashSet<_>>())
    }

    fn gate(name: &str, qubits: &[usize], params: &[f64]) -> Instruction {
        Instruction::gate(name, qubits.to_vec(), params.to_vec())
    }

    fn equivalent(a: &[Instruction], b: &[Instruction], n: usize) -> bool {
        dense_equal_up_to_phase(&circuit_unitary(a, n), &circuit_unitary(b, n), 1e-8)
    }

    #[test]
    fn an_x_travels_through_cx_targets_to_a_run_that_absorbs_it() {
        // As the 8-qubit adder's gathered blocks leave it: an X alone after
        // a CX it is the target of, two CX targets away from a generic run.
        use std::f64::consts::FRAC_PI_4;
        let circuit = vec![
            gate("cx", &[1, 2], &[]),
            gate("x", &[2], &[]),
            gate("rz", &[2], &[FRAC_PI_4]),
            gate("cx", &[0, 2], &[]),
            gate("rz", &[2], &[-FRAC_PI_4]),
            gate("cx", &[1, 2], &[]),
            gate("rz", &[2], &[0.3]),
            gate("sx", &[2], &[]),
            gate("rz", &[2], &[0.7]),
            gate("cx", &[2, 0], &[]),
        ];
        let b = bases(&["rz", "sx", "x", "cx"]);
        let out = insert_pauli_pairs(circuit.clone(), &b);
        assert!(equivalent(&out, &circuit, 3));
        assert_eq!(cost(&circuit), (4, 10));
        assert_eq!(cost(&out), (4, 9));
        assert!(out.iter().all(|i| i.name.as_str() != "x"));
    }

    #[test]
    fn a_pauli_does_not_cross_a_gate_it_does_not_commute_with() {
        // An X on a CX's control would become X on both qubits; a
        // measurement stops every Pauli.
        let b = bases(&["rz", "sx", "x", "cx"]);
        let through_control = vec![
            gate("x", &[0], &[]),
            gate("cx", &[0, 1], &[]),
            gate("rz", &[0], &[0.2]),
            gate("sx", &[0], &[]),
            gate("rz", &[0], &[0.4]),
        ];
        assert_eq!(
            insert_pauli_pairs(through_control.clone(), &b),
            through_control
        );
        let measured = vec![
            gate("x", &[0], &[]),
            Instruction {
                clbits: vec![0],
                ..gate("measure", &[0], &[])
            },
            gate("rz", &[0], &[0.2]),
            gate("sx", &[0], &[]),
            gate("rz", &[0], &[0.4]),
        ];
        assert_eq!(insert_pauli_pairs(measured.clone(), &b), measured);
    }

    #[test]
    fn two_xs_meet_and_cancel_and_a_z_crosses_a_cz() {
        let b = bases(&["rz", "sx", "x", "cx"]);
        let circuit = vec![
            gate("x", &[1], &[]),
            gate("cx", &[0, 1], &[]),
            gate("cx", &[2, 1], &[]),
            gate("x", &[1], &[]),
        ];
        let out = insert_pauli_pairs(circuit.clone(), &b);
        assert!(equivalent(&out, &circuit, 3));
        assert_eq!(cost(&out), (2, 2));
        // With CZ, the Paulis that commute are Z on either qubit: a Z
        // rotation by π joins the next run of Z rotations.
        let b = bases(&["rz", "sx", "x", "cz"]);
        let circuit = vec![
            gate("rz", &[0], &[std::f64::consts::PI]),
            gate("sx", &[0], &[]),
            gate("cz", &[0, 1], &[]),
            gate("rz", &[0], &[std::f64::consts::PI]),
            gate("cz", &[0, 1], &[]),
            gate("rz", &[0], &[0.3]),
            gate("sx", &[0], &[]),
        ];
        let out = insert_pauli_pairs(circuit.clone(), &b);
        assert!(equivalent(&out, &circuit, 2));
        assert!(cost(&out).1 < cost(&circuit).1);
    }

    /// `m` as native gates on `qubit`.
    fn run(m: &M2, qubit: usize, b: &[Basis1q]) -> Vec<Instruction> {
        synthesize_1q(m, b)
            .into_iter()
            .map(|op| Instruction::gate(&op.name, [qubit], op.params))
            .collect()
    }

    #[test]
    fn a_z_rotation_travels_through_cx_controls_and_czs() {
        // The last Z rotation of a generic run on qubit 0, carried through
        // the controls of two CXs and a CZ, merges into the first Z rotation
        // of the next generic run: one gate fewer. No Pauli does that.
        use crate::gates::{rz, u3};
        let b = bases(&["rz", "sx", "x", "cx", "cz"]);
        let (first, second) = (u3(0.7, 0.4, 1.1), u3(1.3, -0.6, 0.2));
        let mut circuit = run(&first, 0, &b);
        circuit.extend([
            gate("cx", &[0, 1], &[]),
            gate("cz", &[0, 2], &[]),
            gate("cx", &[0, 2], &[]),
        ]);
        circuit.extend(run(&second, 0, &b));
        assert_eq!(cost(&circuit), (3, 13));
        let out = insert_rotation_pairs(circuit.clone(), &b);
        assert!(equivalent(&out, &circuit, 3));
        assert_eq!(cost(&out), (3, 12));
        assert_eq!(insert_pauli_pairs(circuit.clone(), &b), circuit);
        // A Z rotation alone between two runs goes into the second, with
        // the first's last one: two gates fewer.
        let mut circuit = run(&first, 0, &b);
        circuit.push(gate("cz", &[0, 1], &[]));
        circuit.extend(run(&rz(0.9), 0, &b));
        circuit.push(gate("cx", &[0, 2], &[]));
        circuit.extend(run(&second, 0, &b));
        let out = insert_rotation_pairs(circuit.clone(), &b);
        assert!(equivalent(&out, &circuit, 3));
        assert_eq!(cost(&out).1, cost(&circuit).1 - 2);
    }

    #[test]
    fn an_x_rotation_travels_through_cx_targets() {
        // Runs on the target of two CXs, a √X each, the first ending and
        // the second starting with X rotations that cancel through both CXs
        // once moved: the pass finds what was put there.
        use crate::gates::{rx, sx};
        let b = bases(&["rz", "sx", "x", "cx"]);
        let (first, second) = (sx(), sx());
        let plain = |a: M2, z: M2| {
            let mut c = run(&a, 2, &b);
            c.extend([gate("cx", &[0, 2], &[]), gate("cx", &[1, 2], &[])]);
            c.extend(run(&z, 2, &b));
            c
        };
        let simple = plain(first, second);
        let hidden = plain(rx(0.8) * first, second * rx(-0.8));
        assert!(equivalent(&simple, &hidden, 3));
        assert!(cost(&hidden).1 > cost(&simple).1);
        let out = insert_rotation_pairs(hidden.clone(), &b);
        assert!(equivalent(&out, &hidden, 3));
        assert!(cost(&out).1 <= cost(&simple).1, "{:?}", cost(&out));
        assert!(cost(&insert_pauli_pairs(hidden, &b)).1 > cost(&simple).1);
    }

    #[test]
    fn rotation_pairs_are_never_costlier_than_pauli_pairs() {
        // Random native circuits, with Z rotations, √X and Paulis between
        // CX-like gates: the rotations keep the unitary and the two-qubit
        // gates, and never cost more than the Paulis alone, which are among
        // their choices; on many they cost less.
        let mut rng = Rng::new(29);
        let mut cheaper = 0;
        for native in [
            ["rz", "sx", "x", "cx"],
            ["rz", "sx", "x", "cz"],
            ["rz", "sx", "x", "ecr"],
            ["rz", "rx", "ry", "cz"],
        ] {
            let mut lowering = Lowering::new(native.iter().map(|s| s.to_string()));
            let b = bases(&native);
            for _ in 0..150 {
                let n = 2 + rng.below(3);
                let mut circuit = Vec::new();
                for _ in 0..20 {
                    let raw = match rng.below(6) {
                        0 => gate("x", &[rng.below(n)], &[]),
                        1 | 2 => gate("rz", &[rng.below(n)], &[rng.range(-3.0, 3.0)]),
                        3 => gate("sx", &[rng.below(n)], &[]),
                        _ => {
                            let a = rng.below(n);
                            gate(native[3], &[a, (a + 1 + rng.below(n - 1)) % n], &[])
                        }
                    };
                    circuit.extend(lowering.lower(&raw).unwrap());
                }
                let paulis = insert_pauli_pairs(circuit.clone(), &b);
                let out = insert_rotation_pairs(circuit.clone(), &b);
                assert!(equivalent(&out, &circuit, n), "{native:?}");
                assert_eq!(cost(&out).0, cost(&circuit).0);
                assert!(cost(&out) <= cost(&circuit));
                assert!(cost(&out).1 <= cost(&paulis).1, "{native:?}");
                if cost(&out).1 < cost(&paulis).1 {
                    cheaper += 1;
                }
            }
        }
        assert!(cheaper > 50, "{cheaper}");
    }

    #[test]
    fn pricing_the_sides_apart_changes_nothing() {
        // The same choices, gate for gate, with and without the sum of
        // `Basis1q::split_count`, on circuits where most runs lie between
        // two gates that let Z rotations pass (CZs, CX controls).
        let mut rng = Rng::new(31);
        let mut changed = 0;
        for native in [["rz", "sx", "x", "cz"], ["rz", "sx", "x", "cx"]] {
            let mut lowering = Lowering::new(native.iter().map(|s| s.to_string()));
            let b = bases(&native);
            for _ in 0..200 {
                let n = 2 + rng.below(3);
                let mut circuit = Vec::new();
                for _ in 0..24 {
                    let q = rng.below(n);
                    let raw = match rng.below(5) {
                        0 => gate("rz", &[q], &[rng.range(-3.0, 3.0)]),
                        1 => gate("sx", &[q], &[]),
                        2 => gate("rz", &[q], &[FRAC_PI_2 * rng.below(4) as f64]),
                        _ => gate(native[3], &[q, (q + 1 + rng.below(n - 1)) % n], &[]),
                    };
                    circuit.extend(lowering.lower(&raw).unwrap());
                }
                let apart = insert_pairs(circuit.clone(), &b, true, true);
                assert_eq!(apart, insert_pairs(circuit.clone(), &b, true, false));
                if apart != circuit {
                    changed += 1;
                }
            }
        }
        assert!(changed > 100, "{changed}");
    }

    #[test]
    fn random_native_circuits_keep_their_unitary_and_never_cost_more() {
        let mut rng = Rng::new(17);
        let mut cheaper = 0;
        for native in [
            ["rz", "sx", "x", "cx"],
            ["rz", "sx", "x", "cz"],
            ["rz", "sx", "x", "ecr"],
        ] {
            let mut lowering = Lowering::new(native.iter().map(|s| s.to_string()));
            let b = bases(&native);
            for _ in 0..150 {
                let n = 2 + rng.below(3);
                let mut circuit = Vec::new();
                for _ in 0..16 {
                    let raw = match rng.below(6) {
                        0 => gate("x", &[rng.below(n)], &[]),
                        1 => gate("rz", &[rng.below(n)], &[rng.range(-3.0, 3.0)]),
                        2 => gate("sx", &[rng.below(n)], &[]),
                        3 => gate("z", &[rng.below(n)], &[]),
                        _ => {
                            let a = rng.below(n);
                            gate(native[3], &[a, (a + 1 + rng.below(n - 1)) % n], &[])
                        }
                    };
                    circuit.extend(lowering.lower(&raw).unwrap());
                }
                let out = insert_pauli_pairs(circuit.clone(), &b);
                assert!(equivalent(&out, &circuit, n), "{native:?}");
                assert!(cost(&out) <= cost(&circuit));
                assert_eq!(cost(&out).0, cost(&circuit).0);
                if cost(&out) < cost(&circuit) {
                    cheaper += 1;
                }
            }
        }
        // Not vacuous: the pairs paid on a good share of the circuits.
        assert!(cheaper > 50, "{cheaper}");
    }
}
