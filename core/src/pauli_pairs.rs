//! Pairs of Paulis around two-qubit gates, chosen wire by wire.
//!
//! If a Pauli `P` on one qubit commutes with a two-qubit gate `G` (X with a
//! CX's target, Z with its control or with either qubit of a CZ), then
//! `G = P G P`: a `P` can be put on that qubit just before `G` and another
//! just after it without changing the circuit. Each `P` merges into the
//! single-qubit run on its side, so the run before `G` becomes `P R` and the
//! one after `R' P`. Doing so around several gates along a wire carries a
//! Pauli from one run to another, however many such gates lie between: an X
//! left behind by a block's synthesis, alone between two CXs it is the
//! target of, can travel to a run further on that absorbs it for free (in
//! {rz, sx, x}, `X·Rz(b)·SX·Rz(a)` is `Rz(π − b)·SX·Rz(π + a)`, still three
//! gates), or meet another X and cancel.
//!
//! Each pair leaves the circuit's unitary as it was, and pairs on different
//! wires do not interact, so every wire is solved on its own and exactly:
//! the wire is a chain of runs separated by the gates on it, the choice at
//! each gate is the Pauli of the pair there (only the identity where none
//! commutes with the gate), and the cost of a run depends only on the
//! choices at the gates on its two sides. A dynamic programme over the
//! chain, four states per gate, finds the choices that make the runs
//! cheapest in all; a wire changes only if that is strictly cheaper than
//! leaving it alone, and only its changed runs are resynthesized.
//!
//! Block resynthesis and the shifts through blocks (`symmetric_blocks.rs`)
//! look at one block and the runs just around it; a Pauli that pays only two
//! or more gates further along the wire is out of their sight. On QCC's
//! final circuits these are mostly the X gates the syntheses leave (the
//! 32-qubit Heisenberg chain had 66).

use crate::gates::{matrix_1q, matrix_2q, pauli};
use crate::hashing::FastMap;
use crate::ir::{Instruction, Name};
use crate::linalg::{M2, M4, equal_up_to_phase4, kron};
use crate::one_qubit::{Basis1q, count_1q, synthesize_1q};

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

/// The masks of [`commuting`] by gate name and parameters: a circuit has a
/// few kinds of two-qubit gate (a CX for most), each tried against three
/// Paulis on each side once.
type Masks = FastMap<(Name, [u64; 3]), [u8; 2]>;

/// For each qubit of `inst`, the Paulis (bit `k` for Pauli `k` of
/// [`paulis`]) that commute with it there, up to a global phase. Only the
/// identity for anything but a two-qubit gate with a known matrix: a
/// measurement, a reset, a barrier or a gate on three qubits stops a
/// Pauli.
fn commuting(inst: &Instruction, paulis: &[M2; 4], known: &mut Masks) -> Vec<u8> {
    if inst.qubits.len() != 2 || !inst.is_gate() {
        return vec![1u8; inst.qubits.len()];
    }
    let mut bits = [0u64; 3];
    for (b, p) in bits.iter_mut().zip(inst.params.as_slice()) {
        *b = p.to_bits();
    }
    let masks = known.entry((inst.name.clone(), bits)).or_insert_with(|| {
        let mut masks = [1u8; 2];
        let Some(g) = matrix_2q(&inst.name, &inst.params) else {
            return masks;
        };
        for (k, mask) in masks.iter_mut().enumerate() {
            for (p, pm) in paulis.iter().enumerate().skip(1) {
                let on: M4 = if k == 0 {
                    kron(pm, &M2::IDENTITY)
                } else {
                    kron(&M2::IDENTITY, pm)
                };
                if equal_up_to_phase4(&(on * g * on), &g, 1e-9) {
                    *mask |= 1 << p;
                }
            }
        }
        masks
    });
    masks.to_vec()
}

/// Insert pairs of Paulis around the two-qubit gates where that makes the
/// single-qubit runs cheaper in `bases` (see the module's description).
/// Changed runs are resynthesized, new ones created where a wire had none;
/// every other run is left as it was.
pub fn insert_pauli_pairs(instructions: Vec<Instruction>, bases: &[Basis1q]) -> Vec<Instruction> {
    if bases.is_empty() {
        return instructions;
    }
    let ps = paulis();
    let mut known = Masks::default();
    // What the product of two Paulis costs where a wire has no run: most
    // slots between two CXs on a wire are empty.
    let alone: [[usize; 4]; 4] =
        std::array::from_fn(|a| std::array::from_fn(|b| count_1q(&(ps[a] * ps[b]), bases)));
    let mut items: Vec<Item> = Vec::with_capacity(instructions.len());
    // Per qubit: its slots so far, the Paulis allowed at the operation that
    // ends each closed slot, and whether its last slot holds a run that is
    // still open.
    let mut slots: Vec<Vec<Slot>> = Vec::new();
    let mut allowed: Vec<Vec<u8>> = Vec::new();
    let grow = |slots: &mut Vec<Vec<Slot>>, allowed: &mut Vec<Vec<u8>>, q: usize| {
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
        for (&q, mask) in inst.qubits.iter().zip(commuting(&inst, &ps, &mut known)) {
            grow(&mut slots, &mut allowed, q);
            slots[q].last_mut().expect("a slot per qubit").until = Some(at);
            slots[q].push(Slot {
                run: None,
                until: None,
            });
            allowed[q].push(mask);
        }
        items.push(Item::Op(inst));
    }

    // New matrices for runs that change, by item index, and runs to create,
    // as (item index to put them before, or None for the end; qubit; matrix).
    let mut changed: Vec<Option<M2>> = vec![None; items.len()];
    let mut created: Vec<(Option<usize>, usize, M2)> = Vec::new();
    for (q, (wire, masks)) in slots.iter().zip(&allowed).enumerate() {
        if masks.iter().all(|&m| m == 1) {
            continue;
        }
        let run_of = |slot: &Slot| match slot.run.map(|r| &items[r]) {
            Some(Item::Run { m, gates, .. }) => (*m, gates.len()),
            _ => (M2::IDENTITY, 0),
        };
        // The cost of slot `s` with Pauli `before` at the gate before it and
        // `after` at the gate after it: as it is when both are the identity.
        let slot_cost = |s: usize, before: usize, after: usize| -> usize {
            if wire[s].run.is_none() {
                return alone[after][before];
            }
            let (m, len) = run_of(&wire[s]);
            if before == 0 && after == 0 {
                len
            } else {
                count_1q(&(ps[after] * m * ps[before]), bases)
            }
        };
        let unchanged: usize = (0..wire.len()).map(|s| slot_cost(s, 0, 0)).sum();
        // best[p]: the cheapest cost of the slots up to the current gate,
        // with Pauli p there; from[g][p]: the Pauli at the gate before.
        let mut best = [usize::MAX; 4];
        best[0] = 0;
        let mut from: Vec<[usize; 4]> = Vec::with_capacity(masks.len());
        for (s, &mask) in masks.iter().enumerate() {
            let mut next = [usize::MAX; 4];
            let mut back = [0usize; 4];
            for (after, slot) in next.iter_mut().enumerate() {
                if mask & (1 << after) == 0 {
                    continue;
                }
                // The identity first and strictly cheaper only: ties keep
                // the earlier choice.
                for (before, &so_far) in best.iter().enumerate() {
                    if so_far == usize::MAX {
                        continue;
                    }
                    let total = so_far + slot_cost(s, before, after);
                    if total < *slot {
                        *slot = total;
                        back[after] = before;
                    }
                }
            }
            best = next;
            from.push(back);
        }
        let last = wire.len() - 1;
        let (total, mut before) = (0..4)
            .filter(|&p| best[p] != usize::MAX)
            .map(|p| (best[p] + slot_cost(last, p, 0), p))
            .min()
            .expect("the identity is always allowed");
        if total >= unchanged {
            continue;
        }
        // Walk back from the end: slot s lies between gate s - 1 and gate s.
        let mut after = 0;
        for s in (0..wire.len()).rev() {
            if before != 0 || after != 0 {
                let (m, _) = run_of(&wire[s]);
                let m = ps[after] * m * ps[before];
                match wire[s].run {
                    Some(r) => changed[r] = Some(m),
                    None => created.push((wire[s].until, q, m)),
                }
            }
            if s > 0 {
                after = before;
                before = from[s - 1][before];
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
