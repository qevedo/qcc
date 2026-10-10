//! Circuit passes: single-qubit merging, commutative cancellation, two-qubit
//! block resynthesis, and the default pipeline. Matches
//! `qevedo/compiler/passes/{optimize,resynthesize,base}.py`.

use crate::gates::{is_known_gate, matrix_1q, matrix_2q};
use crate::ir::{Instruction, NON_GATES, cost};
use crate::linalg::{M2, M4, kron};
use crate::lower::Lowering;
use crate::one_qubit::{Basis1q, remainder, synthesize_1q};
use std::collections::{HashMap, HashSet};
use std::f64::consts::TAU;

// -- single-qubit runs ---------------------------------------------------------------

/// `instructions` with each run of single-qubit gates on a qubit replaced by
/// its shortest native equivalent. A run is kept as it is when it is already
/// native and no shorter. Runs on different qubits commute, so each is emitted
/// where its qubit is next used.
pub fn merge_single_qubit_runs(
    instructions: Vec<Instruction>,
    native: &HashSet<String>,
    bases: &[Basis1q],
) -> Vec<Instruction> {
    let mut result = Vec::with_capacity(instructions.len());
    // Open runs by qubit, with the order they were opened in (for the final flush).
    let mut runs: HashMap<usize, (usize, Vec<Instruction>)> = HashMap::new();
    let mut opened = 0usize;
    let flush = |runs: &mut HashMap<usize, (usize, Vec<Instruction>)>,
                 q: usize,
                 result: &mut Vec<Instruction>| {
        if let Some((_, run)) = runs.remove(&q) {
            result.extend(resynthesize_run(run, bases, native));
        }
    };
    for inst in instructions {
        if inst.qubits.len() == 1
            && !NON_GATES.contains(&inst.name.as_str())
            && is_known_gate(&inst.name)
        {
            let entry = runs.entry(inst.qubits[0]).or_insert_with(|| {
                opened += 1;
                (opened, Vec::new())
            });
            entry.1.push(inst);
            continue;
        }
        for &q in &inst.qubits {
            flush(&mut runs, q, &mut result);
        }
        result.push(inst);
    }
    let mut remaining: Vec<(usize, usize)> =
        runs.iter().map(|(&q, (order, _))| (*order, q)).collect();
    remaining.sort();
    for (_, q) in remaining {
        flush(&mut runs, q, &mut result);
    }
    result
}

fn resynthesize_run(
    run: Vec<Instruction>,
    bases: &[Basis1q],
    native: &HashSet<String>,
) -> Vec<Instruction> {
    if bases.is_empty() {
        return run;
    }
    let mut m = M2::IDENTITY;
    for inst in &run {
        match matrix_1q(&inst.name, &inst.params) {
            Some(g) => m = g * m,
            None => return run,
        }
    }
    let ops = synthesize_1q(&m, bases);
    if run.iter().all(|i| native.contains(&i.name)) && ops.len() >= run.len() {
        return run;
    }
    let q = run[0].qubits[0];
    ops.into_iter()
        .map(|op| Instruction::gate(&op.name, vec![q], op.params))
        .collect()
}

// -- cancellation ----------------------------------------------------------------------

const SELF_INVERSE: [&str; 13] = [
    "x", "y", "z", "h", "id", "cx", "cy", "cz", "ch", "swap", "ecr", "ccx", "cswap",
];
const SYMMETRIC: [&str; 7] = ["cz", "swap", "rxx", "ryy", "rzz", "cp", "cphase"];
const MERGEABLE: [&str; 12] = [
    "rz", "rx", "ry", "p", "u1", "phase", "rzz", "rxx", "ryy", "cp", "cphase", "cu1",
];

fn inverse_of(name: &str) -> Option<&'static str> {
    Some(match name {
        "s" => "sdg",
        "sdg" => "s",
        "t" => "tdg",
        "tdg" => "t",
        "sx" => "sxdg",
        "sxdg" => "sx",
        _ => return None,
    })
}

/// The period of a rotation's angle: a controlled rotation by 2π is a CZ-like phase, not a no-op.
fn rotation_period(name: &str) -> Option<f64> {
    match name {
        "rx" | "ry" | "rz" | "p" | "u1" | "phase" | "rxx" | "ryy" | "rzz" | "rzx" | "cp" => {
            Some(TAU)
        }
        "crx" | "cry" | "crz" => Some(2.0 * TAU),
        _ => None,
    }
}

fn is_inverse_pair(left: &Instruction, right: &Instruction) -> bool {
    let (a, b) = (left.name.as_str(), right.name.as_str());
    if !left.is_gate() || !right.is_gate() {
        return false;
    }
    if left.qubits != right.qubits && !(SYMMETRIC.contains(&a) && a == b) {
        return false;
    }
    if a == b && SELF_INVERSE.contains(&a) {
        return true;
    }
    if inverse_of(a) == Some(b) {
        return true;
    }
    if a == b && left.params.len() == 1 && right.params.len() == 1 {
        if let Some(period) = rotation_period(a) {
            return remainder(left.params[0] + right.params[0], period).abs() < 1e-12;
        }
    }
    false
}

/// How a gate acts on each of its qubits: Z if it commutes with Z there (like a
/// CX control), X if it commutes with X (like a CX target). Two gates commute
/// when every qubit they share has the same letter in both.
fn actions(inst: &Instruction) -> Vec<Option<char>> {
    let n = inst.qubits.len();
    let name = inst.name.as_str();
    let all = |c: char| vec![Some(c); n];
    match name {
        "rz" | "p" | "u1" | "phase" | "z" | "s" | "sdg" | "t" | "tdg" | "id" | "cz" | "cp"
        | "cphase" | "cu1" | "crz" | "rzz" | "ccz" => all('Z'),
        "x" | "sx" | "sxdg" | "rx" | "rxx" => all('X'),
        "cx" | "cnot" | "crx" | "ccx" => {
            let controls = if name == "ccx" { 2 } else { 1 };
            (0..n)
                .map(|i| Some(if i < controls { 'Z' } else { 'X' }))
                .collect()
        }
        "cy" | "ch" | "cry" | "csx" | "cu" | "cu3" | "cswap" => {
            (0..n).map(|i| (i == 0).then_some('Z')).collect()
        }
        _ => vec![None; n],
    }
}

fn commute(a: &Instruction, b: &Instruction) -> bool {
    if !a.is_gate() || !b.is_gate() {
        return !a.qubits.iter().any(|q| b.qubits.contains(q));
    }
    let (actions_a, actions_b) = (actions(a), actions(b));
    for (i, q) in a.qubits.iter().enumerate() {
        if let Some(j) = b.qubits.iter().position(|x| x == q) {
            if actions_a[i].is_none() || actions_a[i] != actions_b[j] {
                return false;
            }
        }
    }
    true
}

enum Combined {
    /// The two gates cancel.
    Cancel,
    /// The two gates are one gate.
    Merged(Instruction),
    NoMerge,
}

/// `left` followed by `right` as one gate.
fn combine(left: &Instruction, right: &Instruction) -> Combined {
    if is_inverse_pair(left, right) {
        return Combined::Cancel;
    }
    let (a, b) = (left.name.as_str(), right.name.as_str());
    if a != b || !MERGEABLE.contains(&a) || left.params.len() != 1 || right.params.len() != 1 {
        return Combined::NoMerge;
    }
    let same_set = left.qubits.len() == right.qubits.len()
        && left.qubits.iter().all(|q| right.qubits.contains(q));
    if left.qubits != right.qubits && !(SYMMETRIC.contains(&a) && same_set) {
        return Combined::NoMerge;
    }
    let angle = left.params[0] + right.params[0];
    let period = rotation_period(a).unwrap_or(TAU);
    if remainder(angle, period).abs() < 1e-12 {
        return Combined::Cancel;
    }
    Combined::Merged(Instruction::gate(
        &left.name,
        left.qubits.clone(),
        vec![angle],
    ))
}

/// Cancel inverse pairs and merge rotations across gates they commute with:
/// `rz(a) q0; cx q0, q1; rz(b) q0` becomes `rz(a + b) q0; cx q0, q1`.
pub fn commutative_cancellation(instructions: Vec<Instruction>) -> Vec<Instruction> {
    /// How far back along a qubit to look for a partner.
    const WINDOW: usize = 32;
    let mut kept: Vec<Option<Instruction>> = Vec::with_capacity(instructions.len());
    let mut timelines: HashMap<usize, Vec<usize>> = HashMap::new();
    for inst in instructions {
        if let Some(partner) = find_partner(&inst, &kept, &timelines, WINDOW) {
            let previous = kept[partner].as_ref().unwrap();
            match combine(previous, &inst) {
                Combined::Cancel => {
                    for q in previous.qubits.clone() {
                        let line = timelines.get_mut(&q).unwrap();
                        let pos = line.iter().position(|&j| j == partner).unwrap();
                        line.remove(pos);
                    }
                    kept[partner] = None;
                    continue;
                }
                Combined::Merged(merged) => {
                    kept[partner] = Some(merged);
                    continue;
                }
                Combined::NoMerge => {}
            }
        }
        for &q in &inst.qubits {
            timelines.entry(q).or_default().push(kept.len());
        }
        kept.push(Some(inst));
    }
    kept.into_iter().flatten().collect()
}

/// An earlier gate on the same qubits that `inst` can reach by commuting.
fn find_partner(
    inst: &Instruction,
    kept: &[Option<Instruction>],
    timelines: &HashMap<usize, Vec<usize>>,
    window: usize,
) -> Option<usize> {
    if inst.qubits.is_empty() || !inst.is_gate() {
        return None;
    }
    let first = timelines.get(&inst.qubits[0])?;
    for &index in first.iter().rev().take(window) {
        let candidate = kept[index].as_ref().unwrap();
        let same_set = candidate.qubits.len() == inst.qubits.len()
            && candidate.qubits.iter().all(|q| inst.qubits.contains(q));
        if same_set
            && !matches!(combine(candidate, inst), Combined::NoMerge)
            && reachable(inst, index, kept, timelines)
        {
            return Some(index);
        }
        if !commute(candidate, inst) {
            return None;
        }
    }
    None
}

/// Whether every gate after `index` on `inst`'s other qubits commutes with `inst`.
fn reachable(
    inst: &Instruction,
    index: usize,
    kept: &[Option<Instruction>],
    timelines: &HashMap<usize, Vec<usize>>,
) -> bool {
    inst.qubits[1..].iter().all(|q| {
        timelines[q]
            .iter()
            .filter(|&&j| j > index)
            .all(|&j| commute(kept[j].as_ref().unwrap(), inst))
    })
}

// -- two-qubit blocks -------------------------------------------------------------------

struct Block {
    pair: [usize; 2],
    instructions: Vec<Instruction>,
}

impl Block {
    /// The block's unitary, with `pair[0]` as the most significant qubit.
    fn matrix(&self) -> Option<M4> {
        let mut out = M4::IDENTITY;
        for inst in &self.instructions {
            let m = if inst.qubits.len() == 1 {
                let g = matrix_1q(&inst.name, &inst.params)?;
                if inst.qubits[0] == self.pair[0] {
                    kron(&g, &M2::IDENTITY)
                } else {
                    kron(&M2::IDENTITY, &g)
                }
            } else {
                let g = matrix_2q(&inst.name, &inst.params)?;
                if inst.qubits[0] == self.pair[0] {
                    g
                } else {
                    g.swapped()
                }
            };
            out = m * out;
        }
        Some(out)
    }
}

enum Item {
    Inst(Instruction),
    Block(usize),
}

/// Replace each maximal run of gates on one pair of qubits by its cheapest
/// synthesis, when that has fewer two-qubit gates, or as many and fewer gates.
pub fn resynthesize_two_qubit_blocks(
    instructions: Vec<Instruction>,
    lowering: &mut Lowering,
) -> Vec<Instruction> {
    // A block stands where its first gate was: gates that come later on other
    // qubits commute with the rest of the block.
    let mut result: Vec<Item> = Vec::new();
    let mut blocks: Vec<Block> = Vec::new();
    let mut open: HashMap<usize, usize> = HashMap::new();
    // Single-qubit gates on qubits in no block yet, with the order they started
    // in; a block that starts on the qubit absorbs them.
    let mut pending: HashMap<usize, (usize, Vec<Instruction>)> = HashMap::new();
    let mut started = 0usize;

    fn close(
        q: usize,
        open: &mut HashMap<usize, usize>,
        blocks: &[Block],
        pending: &mut HashMap<usize, (usize, Vec<Instruction>)>,
        result: &mut Vec<Item>,
    ) {
        if let Some(b) = open.remove(&q) {
            for p in blocks[b].pair {
                open.remove(&p);
            }
        }
        if let Some((_, run)) = pending.remove(&q) {
            result.extend(run.into_iter().map(Item::Inst));
        }
    }

    for inst in instructions {
        let gate = inst.is_gate() && is_known_gate(&inst.name);
        if gate && inst.qubits.len() == 1 {
            let q = inst.qubits[0];
            match open.get(&q) {
                Some(&b) => blocks[b].instructions.push(inst),
                None => {
                    pending
                        .entry(q)
                        .or_insert_with(|| {
                            started += 1;
                            (started, Vec::new())
                        })
                        .1
                        .push(inst);
                }
            }
            continue;
        }
        if gate && inst.qubits.len() == 2 {
            let (a, b) = (inst.qubits[0], inst.qubits[1]);
            if let (Some(&x), Some(&y)) = (open.get(&a), open.get(&b)) {
                if x == y {
                    blocks[x].instructions.push(inst);
                    continue;
                }
            }
            close(a, &mut open, &blocks, &mut pending, &mut result);
            close(b, &mut open, &blocks, &mut pending, &mut result);
            let mut body: Vec<Instruction> = Vec::new();
            for q in [a, b] {
                if let Some((_, run)) = pending.remove(&q) {
                    body.extend(run);
                }
            }
            body.push(inst);
            blocks.push(Block {
                pair: [a, b],
                instructions: body,
            });
            open.insert(a, blocks.len() - 1);
            open.insert(b, blocks.len() - 1);
            result.push(Item::Block(blocks.len() - 1));
            continue;
        }
        for q in inst.qubits.clone() {
            close(q, &mut open, &blocks, &mut pending, &mut result);
        }
        result.push(Item::Inst(inst));
    }
    let mut remaining: Vec<(usize, usize)> =
        pending.iter().map(|(&q, (order, _))| (*order, q)).collect();
    remaining.sort();
    for (_, q) in remaining {
        if let Some((_, run)) = pending.remove(&q) {
            result.extend(run.into_iter().map(Item::Inst));
        }
    }

    let mut blocks: Vec<Option<Block>> = blocks.into_iter().map(Some).collect();
    let mut out = Vec::new();
    for item in result {
        match item {
            Item::Inst(inst) => out.push(inst),
            Item::Block(b) => out.extend(best_for_block(blocks[b].take().unwrap(), lowering)),
        }
    }
    out
}

fn best_for_block(block: Block, lowering: &mut Lowering) -> Vec<Instruction> {
    let two_qubit_gates = block
        .instructions
        .iter()
        .filter(|i| i.qubits.len() == 2)
        .count();
    if two_qubit_gates < 2 && block.instructions.len() <= 3 {
        return block.instructions;
    }
    let Some(m) = block.matrix() else {
        return block.instructions;
    };
    match lowering.two_qubit(&m, block.pair) {
        Ok(candidate) if cost(&candidate) < cost(&block.instructions) => candidate,
        _ => block.instructions,
    }
}

// -- pipeline ---------------------------------------------------------------------------

/// Run named passes in order: "decompose", "merge_1q", "commutative_cancellation"
/// and "resynthesize_2q", sharing one lowering (and its caches).
pub fn run_passes(
    instructions: Vec<Instruction>,
    native: impl IntoIterator<Item = String>,
    passes: &[&str],
) -> Result<Vec<Instruction>, String> {
    let mut lowering = Lowering::new(native);
    let mut current = instructions;
    for pass in passes {
        current = match *pass {
            "decompose" => {
                let mut out = Vec::with_capacity(current.len());
                for inst in &current {
                    out.extend(lowering.lower(inst)?);
                }
                out
            }
            "merge_1q" => merge_single_qubit_runs(current, &lowering.native, &lowering.bases_1q),
            "commutative_cancellation" => commutative_cancellation(current),
            "resynthesize_2q" => resynthesize_two_qubit_blocks(current, &mut lowering),
            other => return Err(format!("unknown pass '{other}'")),
        };
    }
    Ok(current)
}

/// The default pipeline: lower to native gates, then two rounds of
/// [merge single-qubit runs → commutative cancellation → block resynthesis],
/// then a final merge.
pub fn compile(
    instructions: Vec<Instruction>,
    native: impl IntoIterator<Item = String>,
) -> Result<Vec<Instruction>, String> {
    let mut lowering = Lowering::new(native);
    let mut current = Vec::with_capacity(instructions.len());
    for inst in &instructions {
        current.extend(lowering.lower(inst)?);
    }
    for _ in 0..2 {
        current = merge_single_qubit_runs(current, &lowering.native, &lowering.bases_1q);
        current = commutative_cancellation(current);
        current = resynthesize_two_qubit_blocks(current, &mut lowering);
    }
    Ok(merge_single_qubit_runs(
        current,
        &lowering.native,
        &lowering.bases_1q,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gates::signature;
    use crate::testing::{Rng, circuit_unitary, dense_equal_up_to_phase};

    const INPUT: [&str; 33] = [
        "id", "x", "y", "z", "h", "s", "sdg", "t", "tdg", "sx", "sxdg", "rx", "ry", "rz", "p",
        "u1", "u2", "u3", "cx", "cy", "cz", "ch", "csx", "swap", "iswap", "ecr", "cp", "crx",
        "cry", "crz", "rxx", "rzz", "cu",
    ];

    fn random_circuit(rng: &mut Rng, n: usize, length: usize) -> Vec<Instruction> {
        (0..length)
            .map(|_| {
                let name = INPUT[rng.below(INPUT.len())];
                let (arity, count) = signature(name).unwrap();
                let first = rng.below(n);
                let qubits = if arity == 1 {
                    vec![first]
                } else {
                    vec![first, (first + 1 + rng.below(n - 1)) % n]
                };
                Instruction::gate(
                    name,
                    qubits,
                    (0..count).map(|_| rng.range(-3.0, 3.0)).collect(),
                )
            })
            .collect()
    }

    fn natives(names: &[&str]) -> Vec<String> {
        names.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn compiled_circuits_are_equivalent_and_native() {
        let devices: [&[&str]; 4] = [
            &["rz", "sx", "x", "cx"],
            &["rz", "rx", "cz"],
            &["u3", "ecr"],
            &["rx", "ry", "rxx"],
        ];
        let mut rng = Rng::new(9);
        for device in devices {
            for _ in 0..40 {
                let circuit = random_circuit(&mut rng, 3, 25);
                let out = compile(circuit.clone(), natives(device)).unwrap();
                assert!(
                    out.iter().all(|i| device.contains(&i.name.as_str())),
                    "{device:?}"
                );
                assert!(dense_equal_up_to_phase(
                    &circuit_unitary(&out, 3),
                    &circuit_unitary(&circuit, 3),
                    1e-8
                ));
            }
        }
    }

    #[test]
    fn three_qubit_definitions_compile() {
        let circuit = vec![
            Instruction::gate("h", vec![0], vec![]),
            Instruction::gate("ccx", vec![0, 1, 2], vec![]),
            Instruction::gate("cswap", vec![2, 0, 1], vec![]),
        ];
        let out = compile(circuit, natives(&["rz", "sx", "x", "cx"])).unwrap();
        assert!(out.iter().filter(|i| i.name == "cx").count() <= 14);
    }

    #[test]
    fn passes_preserve_the_unitary_and_never_cost_more() {
        let mut rng = Rng::new(4);
        let native = natives(&["rz", "sx", "x", "cx"]);
        let mut lowering = Lowering::new(native.clone());
        for _ in 0..200 {
            let raw = random_circuit(&mut rng, 3, 12);
            let mut circuit = Vec::new();
            for inst in &raw {
                circuit.extend(lowering.lower(inst).unwrap());
            }
            let reference = circuit_unitary(&circuit, 3);
            for out in [
                commutative_cancellation(circuit.clone()),
                resynthesize_two_qubit_blocks(circuit.clone(), &mut lowering),
                merge_single_qubit_runs(circuit.clone(), &lowering.native, &lowering.bases_1q),
            ] {
                assert!(dense_equal_up_to_phase(
                    &circuit_unitary(&out, 3),
                    &reference,
                    1e-8
                ));
                assert!(cost(&out) <= cost(&circuit));
            }
        }
    }

    #[test]
    fn errors() {
        let unknown = vec![Instruction::gate("mystery", vec![0], vec![])];
        assert!(
            compile(unknown, natives(&["rz", "sx", "cx"]))
                .unwrap_err()
                .contains("cannot lower gate 'mystery'")
        );
        let no_2q = vec![Instruction::gate("swap", vec![0, 1], vec![])];
        assert!(
            compile(no_2q, natives(&["rz", "sx", "iswap"]))
                .unwrap_err()
                .contains("no supported two-qubit gate")
        );
    }

    #[test]
    fn commutation_rules() {
        let g =
            |name: &str, q: &[usize], p: &[f64]| Instruction::gate(name, q.to_vec(), p.to_vec());
        let out = commutative_cancellation(vec![
            g("rz", &[0], &[0.25]),
            g("cx", &[0, 1], &[]),
            g("rz", &[0], &[0.5]),
        ]);
        assert_eq!(out, vec![g("rz", &[0], &[0.75]), g("cx", &[0, 1], &[])]);
        let blocked = commutative_cancellation(vec![
            g("cx", &[0, 1], &[]),
            g("cx", &[1, 2], &[]),
            g("cx", &[0, 1], &[]),
        ]);
        assert_eq!(blocked.len(), 3);
        let twice = commutative_cancellation(vec![
            g("crx", &[0, 1], &[std::f64::consts::PI]),
            g("crx", &[0, 1], &[std::f64::consts::PI]),
        ]);
        // crx(π) crx(π) = crx(2π) is a CZ-like phase, not the identity: both stay.
        assert_eq!(twice.len(), 2);
    }
}
