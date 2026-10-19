//! Moving single-qubit gates through symmetric two-qubit blocks.
//!
//! A block of gates on one pair of qubits is `U = (A1 ⊗ A0) C (B1 ⊗ B0)`, with
//! `C` canonical (see `two_qubit.rs`). When the Weyl point of `C` has
//! a = b = |c|, as for XX + YY + ZZ, `C` commutes with `V1 ⊗ V0` for every
//! single-qubit `V`, where `V1 = V` and `V0 = V`, or `Z V Z` when c < 0
//! (XX + YY − ZZ is XX + YY + ZZ conjugated by a Z on one qubit). So
//! `U X = Y U` with `X_q = B_q† V_q B_q` and `Y_q = A_q V_q A_q†`. The single-qubit runs around the block can then trade
//! a factor through it: `R_after U R_before = (R_after Y) U (X† R_before)`.
//!
//! Each block is synthesized without seeing its neighbours' runs, so the
//! runs between symmetric blocks are generic. This pass picks `V` for each
//! such block, in order, to make the four runs around it cheapest: on a
//! Heisenberg chain compiled for CZ, the runs between blocks are what made
//! QCC's output a fifth longer than Qiskit's.

use crate::gates::{matrix_1q, pauli};
use crate::ir::Instruction;
use crate::linalg::{M2, equal_up_to_phase4, kron};
use crate::one_qubit::{Basis1q, count_1q, synthesize_1q};
use crate::passes::Block;
use crate::two_qubit::{TOL, kak};

/// A run of single-qubit gates on one qubit, or any other instruction.
enum Item {
    Run {
        qubit: usize,
        m: M2,
        gates: Vec<Instruction>,
        changed: bool,
    },
    Op(Instruction),
}

/// A maximal run of two-qubit gates on one pair, by item index.
struct SymBlock {
    pair: [usize; 2],
    ops: Vec<usize>,
    /// The run just before the first gate, and just after the last, on each
    /// qubit of `pair`, if there is one.
    before: [Option<usize>; 2],
    after: [Option<usize>; 2],
}

/// Move single-qubit factors through the blocks whose two-qubit core is
/// symmetric under `V ⊗ V` (or `V ⊗ Z V Z`), where that makes the runs around them cheaper in
/// `bases`. Runs that change are resynthesized; the others are left alone.
pub fn shift_through_symmetric_blocks(
    instructions: Vec<Instruction>,
    bases: &[Basis1q],
) -> Vec<Instruction> {
    if bases.is_empty() {
        return instructions;
    }
    let mut items = into_items(instructions);
    let blocks = find_blocks(&items);
    // Runs to create, as (item index to put them before, qubit, matrix).
    let mut created: Vec<(usize, usize, M2)> = Vec::new();
    for block in &blocks {
        shift(&mut items, block, bases, &mut created);
    }
    created.sort_by_key(|&(at, _, _)| at);
    let mut created = created.into_iter().peekable();
    let mut out = Vec::new();
    for (i, item) in items.into_iter().enumerate() {
        while let Some((_, qubit, m)) = created.next_if(|&(at, _, _)| at == i) {
            out.extend(emit(qubit, &m, bases));
        }
        match item {
            Item::Run {
                qubit,
                m,
                changed: true,
                ..
            } => out.extend(emit(qubit, &m, bases)),
            Item::Run { gates, .. } => out.extend(gates),
            Item::Op(inst) => out.push(inst),
        }
    }
    for (_, qubit, m) in created {
        out.extend(emit(qubit, &m, bases));
    }
    out
}

fn emit(qubit: usize, m: &M2, bases: &[Basis1q]) -> Vec<Instruction> {
    synthesize_1q(m, bases)
        .into_iter()
        .map(|op| Instruction::gate(&op.name, vec![qubit], op.params))
        .collect()
}

/// Group consecutive single-qubit gates on each qubit into runs, each placed
/// where its first gate was.
fn into_items(instructions: Vec<Instruction>) -> Vec<Item> {
    let mut items = Vec::new();
    let mut open: Vec<Option<usize>> = Vec::new();
    for inst in instructions {
        for &q in &inst.qubits {
            if open.len() <= q {
                open.resize(q + 1, None);
            }
        }
        let matrix = match inst.qubits.as_slice() {
            [_] if inst.is_gate() => matrix_1q(&inst.name, &inst.params),
            _ => None,
        };
        match matrix {
            Some(g) => {
                let q = inst.qubits[0];
                match open[q] {
                    Some(r) => {
                        if let Item::Run { m, gates, .. } = &mut items[r] {
                            *m = g * *m;
                            gates.push(inst);
                        }
                    }
                    None => {
                        open[q] = Some(items.len());
                        items.push(Item::Run {
                            qubit: q,
                            m: g,
                            gates: vec![inst],
                            changed: false,
                        });
                    }
                }
            }
            None => {
                for &q in &inst.qubits {
                    open[q] = None;
                }
                items.push(Item::Op(inst));
            }
        }
    }
    items
}

/// The maximal runs of two-qubit gates on one pair, with the runs around them.
fn find_blocks(items: &[Item]) -> Vec<SymBlock> {
    let mut blocks: Vec<SymBlock> = Vec::new();
    // Per qubit: the block it is in, the run since its last operation, and
    // the block whose gate was its last operation.
    let mut open: Vec<Option<usize>> = Vec::new();
    let mut run: Vec<Option<usize>> = Vec::new();
    let mut last: Vec<Option<usize>> = Vec::new();
    let grow = |v: &mut Vec<Option<usize>>, q: usize| {
        if v.len() <= q {
            v.resize(q + 1, None);
        }
    };
    for (i, item) in items.iter().enumerate() {
        match item {
            Item::Run { qubit, .. } => {
                let q = *qubit;
                grow(&mut run, q);
                grow(&mut last, q);
                run[q] = Some(i);
                if let Some(b) = last[q] {
                    let k = blocks[b].pair.iter().position(|&p| p == q).unwrap();
                    blocks[b].after[k] = Some(i);
                }
            }
            Item::Op(inst) => {
                for &q in &inst.qubits {
                    grow(&mut open, q);
                    grow(&mut run, q);
                    grow(&mut last, q);
                }
                let two_qubit = inst.qubits.len() == 2 && inst.is_gate();
                let (a, b) = if two_qubit {
                    (inst.qubits[0], inst.qubits[1])
                } else {
                    (usize::MAX, usize::MAX)
                };
                let same = two_qubit && open[a].is_some() && open[a] == open[b];
                if same {
                    let k = open[a].unwrap();
                    blocks[k].ops.push(i);
                    // The runs since the last gate are inside the block.
                    blocks[k].after = [None, None];
                } else if two_qubit {
                    for q in [a, b] {
                        if let Some(k) = open[q] {
                            for p in blocks[k].pair {
                                open[p] = None;
                            }
                        }
                    }
                    blocks.push(SymBlock {
                        pair: [a, b],
                        ops: vec![i],
                        before: [run[a], run[b]],
                        after: [None, None],
                    });
                    open[a] = Some(blocks.len() - 1);
                    open[b] = Some(blocks.len() - 1);
                } else {
                    for &q in &inst.qubits {
                        if let Some(k) = open[q] {
                            for p in blocks[k].pair {
                                open[p] = None;
                            }
                        }
                    }
                }
                for &q in &inst.qubits {
                    run[q] = None;
                    last[q] = if two_qubit { open[q] } else { None };
                }
            }
        }
    }
    blocks
}

fn run_matrix(items: &[Item], index: Option<usize>) -> M2 {
    match index.map(|i| &items[i]) {
        Some(Item::Run { m, .. }) => *m,
        _ => M2::IDENTITY,
    }
}

/// Choose `V` for one block and update the runs around it.
fn shift(
    items: &mut [Item],
    block: &SymBlock,
    bases: &[Basis1q],
    created: &mut Vec<(usize, usize, M2)>,
) {
    // The block's core: its gates and the runs between them.
    let (first, last) = (block.ops[0], *block.ops.last().unwrap());
    let mut core = Vec::new();
    for item in &items[first..=last] {
        match item {
            Item::Op(inst) => core.push(inst.clone()),
            Item::Run { qubit, gates, .. } if block.pair.contains(qubit) => {
                core.extend(gates.iter().cloned())
            }
            Item::Run { .. } => {}
        }
    }
    let Some(u) = (Block {
        pair: block.pair,
        instructions: core,
    })
    .matrix() else {
        return;
    };
    let k = kak(&u);
    let [a, b, c] = k.point;
    if a < TOL || (a - b).abs() > TOL || (b - c.abs()).abs() > TOL {
        return;
    }
    // V on each qubit of the pair, from V on the first.
    let z = pauli(2);
    let on = |v: &M2| -> [M2; 2] { [*v, if c < 0.0 { z * *v * z } else { *v }] };
    // Qubit pair[0] is the first (most significant) factor: a1 and b1.
    let (outer_a, outer_b) = ([k.a1, k.a0], [k.b1, k.b0]);
    let before = block.before.map(|r| run_matrix(items, r));
    let after = block.after.map(|r| run_matrix(items, r));
    let factors = |v: &M2| -> ([M2; 2], [M2; 2]) {
        let v = on(v);
        (
            [0, 1].map(|j| outer_b[j].adjoint() * v[j] * outer_b[j]),
            [0, 1].map(|j| outer_a[j] * v[j] * outer_a[j].adjoint()),
        )
    };
    let moved = |v: &M2| -> [M2; 4] {
        let (x, y) = factors(v);
        [
            x[0].adjoint() * before[0],
            x[1].adjoint() * before[1],
            after[0] * y[0],
            after[1] * y[1],
        ]
    };
    let total = |runs: &[M2; 4]| -> usize { runs.iter().map(|m| count_1q(m, bases)).sum() };
    let current = total(&[before[0], before[1], after[0], after[1]]);
    // V that turns one of the four runs into the identity (`on` is its own
    // inverse on each qubit).
    let mut candidates = Vec::new();
    for j in 0..2 {
        for v in [
            outer_b[j] * before[j] * outer_b[j].adjoint(),
            outer_a[j].adjoint() * after[j].adjoint() * outer_a[j],
        ] {
            candidates.push(on(&v)[j]);
        }
    }
    let Some((best, v)) = candidates
        .into_iter()
        .map(|v| (total(&moved(&v)), v))
        .min_by_key(|(t, _)| *t)
    else {
        return;
    };
    if best >= current {
        return;
    }
    // Y U X† = U, checked rather than trusted to the tolerance above.
    let (x, y) = factors(&v);
    let shifted = kron(&y[0], &y[1]) * u * kron(&x[0], &x[1]).adjoint();
    if !equal_up_to_phase4(&shifted, &u, 1e-8) {
        return;
    }
    let runs = moved(&v);
    let slots = [
        (block.before[0], first, block.pair[0]),
        (block.before[1], first, block.pair[1]),
        (block.after[0], last + 1, block.pair[0]),
        (block.after[1], last + 1, block.pair[1]),
    ];
    for ((index, at, qubit), m) in slots.into_iter().zip(runs) {
        match index {
            Some(r) => {
                if let Item::Run {
                    m: old, changed, ..
                } = &mut items[r]
                {
                    *old = m;
                    *changed = true;
                }
            }
            None if !m.close(&M2::IDENTITY, 1e-9) => created.push((at, qubit, m)),
            None => {}
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::passes::compile;
    use crate::testing::{Rng, circuit_unitary, dense_equal_up_to_phase};

    fn heisenberg(n: usize, steps: usize, rng: &mut Rng) -> Vec<Instruction> {
        let mut out = Vec::new();
        for _ in 0..steps {
            for parity in [0, 1] {
                for q in (parity..n - 1).step_by(2) {
                    let t = rng.range(0.05, 0.5);
                    for name in ["rxx", "ryy", "rzz"] {
                        out.push(Instruction::gate(name, vec![q, q + 1], vec![t]));
                    }
                }
            }
            // Fields between the steps, so that the runs are not trivial.
            for q in 0..n {
                out.push(Instruction::gate("rx", vec![q], vec![rng.range(-1.0, 1.0)]));
            }
        }
        out
    }

    #[test]
    fn shifting_keeps_the_unitary_and_saves_gates() {
        let mut rng = Rng::new(3);
        let native: Vec<String> = ["rz", "sx", "x", "cz"].map(String::from).to_vec();
        for _ in 0..10 {
            let circuit = heisenberg(4, 3, &mut rng);
            let compiled = compile(circuit.clone(), native.clone()).unwrap();
            assert!(dense_equal_up_to_phase(
                &circuit_unitary(&compiled, 4),
                &circuit_unitary(&circuit, 4),
                1e-7
            ));
        }
    }

    #[test]
    fn a_symmetric_block_takes_its_neighbours_gates() {
        // Rz Rx on both qubits, then XX + YY + ZZ, then their inverses: with
        // V moved through the block, the outer runs cancel to nothing.
        let native: std::collections::HashSet<String> =
            ["rz", "sx", "x", "cx"].map(String::from).into();
        let bases = crate::one_qubit::bases_for(&native);
        let merged = |c: Vec<Instruction>| {
            crate::passes::merge_single_qubit_runs(c, &native, &bases)
                .iter()
                .filter(|i| i.qubits.len() == 1)
                .count()
        };
        let g = |name: &str, q: Vec<usize>, p: f64| Instruction::gate(name, q, vec![p]);
        let mut circuit = Vec::new();
        for q in [0, 1] {
            circuit.push(g("rz", vec![q], 0.7));
            circuit.push(g("rx", vec![q], 0.4));
        }
        for name in ["rxx", "ryy", "rzz"] {
            circuit.push(g(name, vec![0, 1], 0.3));
        }
        for q in [0, 1] {
            circuit.push(g("rx", vec![q], -0.4));
            circuit.push(g("rz", vec![q], -0.7));
        }
        let out = shift_through_symmetric_blocks(circuit.clone(), &bases);
        assert!(merged(out.clone()) < merged(circuit.clone()));
        assert!(dense_equal_up_to_phase(
            &circuit_unitary(&out, 2),
            &circuit_unitary(&circuit, 2),
            1e-8
        ));
    }

    #[test]
    fn a_lone_heisenberg_block_has_cheap_outer_gates() {
        // The KAK decomposition's outer gates are arbitrary at a = b = |c|;
        // two_qubit.rs tries the symmetric ones. Qiskit's level 3 uses 20.
        let circuit: Vec<Instruction> = ["rxx", "ryy", "rzz"]
            .iter()
            .map(|name| Instruction::gate(name, vec![0, 1], vec![0.2]))
            .collect();
        let native: Vec<String> = ["rz", "sx", "x", "cz"].map(String::from).to_vec();
        let out = compile(circuit.clone(), native).unwrap();
        assert!(out.len() <= 20, "{} gates", out.len());
        assert!(dense_equal_up_to_phase(
            &circuit_unitary(&out, 2),
            &circuit_unitary(&circuit, 2),
            1e-8
        ));
    }
}
