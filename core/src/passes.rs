//! Circuit passes: single-qubit merging, commutative cancellation, two-qubit
//! block resynthesis, and the default pipeline (which also runs
//! `phase_folding`). Matches
//! `qevedo/compiler/passes/{optimize,resynthesize,base}.py`.

use crate::gates::{is_known_gate, matrix_1q, matrix_2q, signature};
use crate::ir::{Instruction, NON_GATES, cost};
use crate::linalg::{M2, M4, kron};
use crate::lower::Lowering;
use crate::one_qubit::{Basis1q, count_1q, remainder, synthesize_1q};
use crate::permutation::measured_out;
use crate::phase_folding::phase_folding;
use crate::symmetric_blocks::shift_through_symmetric_blocks;
use crate::two_qubit::{Objective, cx_count};
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

/// `cost(&merge_single_qubit_runs(instructions, native, bases))`, without
/// building the merged circuit.
fn merged_cost(
    instructions: &[Instruction],
    native: &HashSet<String>,
    bases: &[Basis1q],
) -> (usize, usize) {
    let mut two_qubit = 0;
    let mut total = 0;
    // Open runs by qubit: their unitary (None once a gate has no matrix),
    // length, and whether they are all native.
    let mut runs: HashMap<usize, (Option<M2>, usize, bool)> = HashMap::new();
    let close = |run: (Option<M2>, usize, bool)| -> usize {
        match run {
            (Some(m), len, all_native) if !bases.is_empty() => {
                let count = count_1q(&m, bases);
                if all_native { count.min(len) } else { count }
            }
            (_, len, _) => len,
        }
    };
    for inst in instructions {
        if inst.qubits.len() == 1
            && !NON_GATES.contains(&inst.name.as_str())
            && is_known_gate(&inst.name)
        {
            let run = runs
                .entry(inst.qubits[0])
                .or_insert((Some(M2::IDENTITY), 0, true));
            run.0 = match (run.0, matrix_1q(&inst.name, &inst.params)) {
                (Some(m), Some(g)) => Some(g * m),
                _ => None,
            };
            run.1 += 1;
            run.2 &= native.contains(&inst.name);
            continue;
        }
        for q in &inst.qubits {
            if let Some(run) = runs.remove(q) {
                total += close(run);
            }
        }
        total += 1;
        if inst.qubits.len() > 1 {
            two_qubit += 1;
        }
    }
    total += runs.into_values().map(close).sum::<usize>();
    (two_qubit, total)
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

#[derive(Clone)]
pub(crate) struct Block {
    pub(crate) pair: [usize; 2],
    pub(crate) instructions: Vec<Instruction>,
}

impl Block {
    /// The block's unitary, with `pair[0]` as the most significant qubit.
    pub(crate) fn matrix(&self) -> Option<M4> {
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
/// synthesis, when that has fewer two-qubit gates, or as many and fewer
/// gates inside the block ([`block_cost`]). As counting only inside can, in
/// rare circuits, cost more once single-qubit runs merge, the result is
/// compared, merged, with replacing only blocks that are cheaper in all
/// their gates and with the input, and the cheapest is kept.
pub fn resynthesize_two_qubit_blocks(
    instructions: Vec<Instruction>,
    lowering: &mut Lowering,
) -> Vec<Instruction> {
    resynthesize_permuting(instructions, lowering, false).0
}

/// [`resynthesize_two_qubit_blocks`], and with `permute`, also permutation-aware:
/// a block U may be synthesized as SWAP·U followed by a SWAP that is elided
/// by relabelling the qubits after it, when that is cheaper (a CX each way is
/// a SWAP and one CX). Returns the instructions and where each qubit's state
/// ends: on wire `wires[q]`.
pub fn resynthesize_permuting(
    instructions: Vec<Instruction>,
    lowering: &mut Lowering,
    permute: bool,
) -> (Vec<Instruction>, Vec<usize>) {
    let qubits = qubit_count(&instructions);
    let (items, blocks) = collect_blocks(instructions.clone());
    // Both measures' syntheses at once, on all cores.
    lowering.prefetch_two_qubit(&synthesis_requests(&blocks, permute));
    let inside = resynthesize_blocks(&items, &blocks, qubits, lowering, Measure::Inside, permute);
    let whole = resynthesize_blocks(&items, &blocks, qubits, lowering, Measure::Whole, permute);
    let unchanged = (0..qubits).collect();
    let merged = |c: &[Instruction]| merged_cost(c, &lowering.native, &lowering.bases_1q);
    // The first of the cheapest, once merged: never worse than the input.
    [inside, whole, (instructions, unchanged)]
        .into_iter()
        .min_by_key(|(c, _)| merged(c))
        .expect("three candidates")
}

/// One more than the highest qubit index of `instructions`.
fn qubit_count(instructions: &[Instruction]) -> usize {
    instructions
        .iter()
        .flat_map(|i| i.qubits.iter().map(|&q| q + 1))
        .max()
        .unwrap_or(0)
}

/// How [`resynthesize_blocks`] compares a block with its synthesis.
#[derive(Clone, Copy)]
enum Measure {
    /// By [`block_cost`]: two-qubit gates, then gates inside the block.
    Inside,
    /// By [`cost`]: two-qubit gates, then all gates.
    Whole,
}

/// The circuit as instructions and maximal blocks of gates on one pair of
/// qubits (with the single-qubit gates before them), in order.
fn collect_blocks(instructions: Vec<Instruction>) -> (Vec<Item>, Vec<Block>) {
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
    (result, blocks)
}

/// Each block replaced by the cheapest of it and its syntheses, by `measure`.
fn resynthesize_blocks(
    items: &[Item],
    blocks: &[Block],
    qubits: usize,
    lowering: &mut Lowering,
    measure: Measure,
    permute: bool,
) -> (Vec<Instruction>, Vec<usize>) {
    // A block synthesized with its SWAP elided moves its qubits' states to
    // each other's wires, and everything after is relabelled to follow them.
    let mut wires: Vec<usize> = (0..qubits).collect();
    let mut out = Vec::new();
    for item in items {
        match item {
            Item::Inst(inst) => out.push(inst.relabel(&wires)),
            &Item::Block(b) => {
                let block = blocks[b].clone();
                let [p, q] = block.pair;
                let (synthesized, swapped) = best_for_block(block, lowering, measure, permute);
                out.extend(synthesized.iter().map(|i| i.relabel(&wires)));
                if swapped {
                    wires.swap(p, q);
                }
            }
        }
    }
    (out, wires)
}

/// The unitaries [`resynthesize_blocks`] will synthesize, under both
/// measures' objectives.
fn synthesis_requests(blocks: &[Block], permute: bool) -> Vec<(M4, Objective)> {
    let swap = matrix_2q("swap", &[]).expect("a standard gate");
    let mut requests = Vec::new();
    for block in blocks {
        if !worth_resynthesizing(block) {
            continue;
        }
        if let Some(m) = block.matrix() {
            let swapped = permute && cx_count(&(swap * m)) < cx_count(&m);
            for objective in [Objective::InCircuit, Objective::Alone] {
                requests.push((m, objective));
                if swapped {
                    requests.push((swap * m, objective));
                }
            }
        }
    }
    requests
}

/// Whether a synthesis could beat the block: not with fewer than two
/// two-qubit gates and three gates in all.
fn worth_resynthesizing(block: &Block) -> bool {
    let two_qubit_gates = block
        .instructions
        .iter()
        .filter(|i| i.qubits.len() == 2)
        .count();
    two_qubit_gates >= 2 || block.instructions.len() > 3
}

/// The cheapest of the block, its synthesis and, with `permute`, the
/// synthesis of SWAP·U; and whether that last one, whose SWAP the caller
/// elides, was chosen.
fn best_for_block(
    block: Block,
    lowering: &mut Lowering,
    measure: Measure,
    permute: bool,
) -> (Vec<Instruction>, bool) {
    if !worth_resynthesizing(&block) {
        return (block.instructions, false);
    }
    let Some(m) = block.matrix() else {
        return (block.instructions, false);
    };
    let better = |candidate: &[Instruction], original: &[Instruction]| match measure {
        Measure::Inside => block_cost(candidate) < block_cost(original),
        Measure::Whole => cost(candidate) < cost(original),
    };
    // Inside counts the gates between the block's first and last two-qubit
    // gates first, and so should the synthesis.
    let objective = match measure {
        Measure::Inside => Objective::InCircuit,
        Measure::Whole => Objective::Alone,
    };
    let mut best = match lowering.two_qubit(&m, block.pair, objective) {
        Ok(candidate) if better(&candidate, &block.instructions) => candidate,
        _ => block.instructions,
    };
    let mut swapped = false;
    // Only when SWAP·U needs fewer CX-like gates than U: otherwise the
    // relabelling saves nothing, and synthesizing every generic block twice
    // would double the time.
    let swap = matrix_2q("swap", &[]).expect("a standard gate");
    if permute && cx_count(&(swap * m)) < cx_count(&m) {
        if let Ok(candidate) = lowering.two_qubit(&(swap * m), block.pair, objective) {
            if better(&candidate, &best) {
                (best, swapped) = (candidate, true);
            }
        }
    }
    (best, swapped)
}

/// What a block costs once its outer single-qubit gates have merged with
/// their neighbours' (as the next single-qubit merge does): its two-qubit
/// gates, then the single-qubit gates between its first and last two-qubit
/// gate, then all its gates. Counting the outer gates in full would refuse a
/// synthesis that saves many gates inside the block for a few more outside.
fn block_cost(instructions: &[Instruction]) -> (usize, usize, usize) {
    let two_qubit: Vec<usize> = (0..instructions.len())
        .filter(|&i| instructions[i].qubits.len() == 2)
        .collect();
    let inside = match (two_qubit.first(), two_qubit.last()) {
        (Some(&first), Some(&last)) => instructions[first..last]
            .iter()
            .filter(|i| i.qubits.len() == 1)
            .count(),
        _ => 0,
    };
    (two_qubit.len(), inside, instructions.len())
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
            "phase_folding" => phase_folding(current),
            "merge_rotations" => crate::pauli_network::merge_rotations(current),
            "pauli_network" => match crate::pauli_network::synthesize(&current, false) {
                Some((out, _)) => out,
                None => current,
            },
            "pauli_network_permuting" => match crate::pauli_network::synthesize(&current, true) {
                Some((out, _)) => out,
                None => current,
            },
            "symmetric_blocks" => shift_through_symmetric_blocks(current, &lowering.bases_1q),
            other => return Err(format!("unknown pass '{other}'")),
        };
    }
    Ok(current)
}

/// The default pipeline: lower to native gates; phase folding, merging
/// single-qubit runs, commutative cancellation and block resynthesis; the
/// same without resynthesis; moving gates through symmetric blocks; and a
/// final merge.
///
/// Every two-qubit gate of the result acts on two qubits that some gate of
/// the input acts on together, so that an input already routed for a device
/// stays routed: this is what the Qiskit plugin needs.
pub fn compile(
    instructions: Vec<Instruction>,
    native: impl IntoIterator<Item = String>,
) -> Result<Vec<Instruction>, String> {
    let native: Vec<String> = native.into_iter().collect();
    Ok(best_pipeline(&instructions, &native, Relabel::Never, false)?.0)
}

/// Whether the pipeline may leave qubits' states on other wires at the end.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Relabel {
    /// No: the output implements the input's unitary.
    Never,
    /// Only when that cannot be seen: every moved state is measured last.
    IfMeasured,
    /// Yes; the caller is told where each state ends.
    Always,
}

/// [`compile`], with permutation-aware block resynthesis when `relabel`
/// allows it. Returns the instructions and where each qubit's state ends: on
/// wire `wires[q]`.
///
/// For a device whose two-qubit gate is CX-like but not CX (CZ, ECR, …), the
/// circuit is also compiled for CX first and then for the device, and the
/// cheaper result kept: phase folding sees parities through CNOTs, but not
/// through the single-qubit gates that turn a CZ into one.
pub fn compile_relabelling(
    instructions: Vec<Instruction>,
    native: impl IntoIterator<Item = String>,
    relabel: Relabel,
) -> Result<(Vec<Instruction>, Vec<usize>), String> {
    let native: Vec<String> = native.into_iter().collect();
    best_pipeline(&instructions, &native, relabel, true)
}

/// The cheapest of [`pipeline`], for CX-like devices [`pipeline`] through
/// CX, and with `new_pairs` (two-qubit gates may act on pairs of qubits the
/// input does not use) the Pauli network, among those `relabel` allows.
pub(crate) fn best_pipeline(
    instructions: &[Instruction],
    native: &[String],
    relabel: Relabel,
    new_pairs: bool,
) -> Result<(Vec<Instruction>, Vec<usize>), String> {
    let permute = relabel != Relabel::Never;
    let mut candidates = vec![pipeline(instructions, native, permute)?];
    if let Some(via) = through_cx(native) {
        if let Ok((middle, first)) = pipeline(instructions, &via, permute) {
            let (out, second) = pipeline(&middle, native, permute)?;
            let wires = first
                .iter()
                .map(|&w| second.get(w).copied().unwrap_or(w))
                .collect();
            candidates.push((out, wires));
        }
    }
    if new_pairs {
        // The network can trade single-qubit gates for CXs, and on the
        // 32-qubit Hamiltonian it would: 663 CXs and 1598 gates in all,
        // against 740 and 1568. Whether that pays depends on the device's
        // error rates, which the pipeline does not know, so it keeps the
        // network only when it needs no more of either.
        let best = candidates.iter().map(|(out, _)| cost(out)).min();
        candidates.extend(through_pauli_network(instructions, native, permute).filter(
            |(out, _)| {
                let (two, all) = cost(out);
                best.is_none_or(|(b_two, b_all)| two <= b_two && all <= b_all)
            },
        ));
    }
    let allowed = |(out, wires): &(Vec<Instruction>, Vec<usize>)| {
        let moved = wires.iter().enumerate().any(|(q, &w)| q != w);
        !moved || relabel == Relabel::Always || measured_out(out, wires)
    };
    // The first of the cheapest: the direct compilation on ties.
    match candidates
        .into_iter()
        .filter(allowed)
        .min_by_key(|(out, _)| cost(out))
    {
        Some(best) => Ok(best),
        None => best_pipeline(instructions, native, Relabel::Never, new_pairs),
    }
}

/// The circuit lowered to CX, synthesized again as a Pauli network
/// (`pauli_network.rs`), then through the [`pipeline`]; None when the
/// circuit is not Clifford gates and Z rotations once lowered, or has too
/// many rotations.
fn through_pauli_network(
    instructions: &[Instruction],
    native: &[String],
    permute: bool,
) -> Option<(Vec<Instruction>, Vec<usize>)> {
    let lowered = lowered_for_network(instructions, native)?;
    let trial = |seed: u64| -> Option<(Vec<Instruction>, Vec<usize>)> {
        let (network, first) = crate::pauli_network::synthesize_seeded(&lowered, permute, seed)?;
        let (out, second) = pipeline(&network, native, permute).ok()?;
        let wires = first
            .iter()
            .map(|&w| second.get(w).copied().unwrap_or(w))
            .collect();
        Some((out, wires))
    };
    let seeds = 0..NETWORK_TRIALS;
    let results: Vec<Option<(Vec<Instruction>, Vec<usize>)>> = if cfg!(target_arch = "wasm32") {
        seeds.map(trial).collect()
    } else {
        // Trials are independent: one thread each.
        std::thread::scope(|scope| {
            let handles: Vec<_> = seeds.map(|seed| scope.spawn(move || trial(seed))).collect();
            handles
                .into_iter()
                .map(|h| h.join().expect("a Pauli-network thread"))
                .collect()
        })
    };
    // The first of the cheapest, so that the result does not depend on the
    // threads.
    results
        .into_iter()
        .flatten()
        .min_by_key(|(out, _)| cost(out))
}

/// How many Pauli networks [`through_pauli_network`] synthesizes, with seeds
/// 0 (no randomness) to this minus one, side by side.
const NETWORK_TRIALS: u64 = 8;

/// The circuit lowered to CX, if the Pauli network is worth trying on it.
fn lowered_for_network(
    instructions: &[Instruction],
    native: &[String],
) -> Option<Vec<Instruction>> {
    let via = if native.iter().any(|g| g == "cx") {
        native.to_vec()
    } else {
        through_cx(native)?
    };
    let mut lowering = Lowering::new(via);
    let mut lowered = Vec::with_capacity(instructions.len());
    for inst in instructions {
        lowered.extend(lowering.lower(inst).ok()?);
    }
    crate::pauli_network::promising(&lowered).then_some(lowered)
}

/// Whether [`best_pipeline`] would try a Pauli network on the circuit.
pub(crate) fn network_promising(instructions: &[Instruction], native: &[String]) -> bool {
    lowered_for_network(instructions, native).is_some()
}

/// `native` with CX as its two-qubit gate, if its two-qubit gates are
/// CX-like but do not include CX.
fn through_cx(native: &[String]) -> Option<Vec<String>> {
    let two_qubit = |name: &str| signature(name).is_some_and(|(arity, _)| arity == 2);
    let cx_like = ["cz", "cy", "ch", "ecr"];
    let gates: Vec<&String> = native.iter().filter(|g| two_qubit(g)).collect();
    if gates.is_empty()
        || native.iter().any(|g| g == "cx")
        || !gates.iter().all(|g| cx_like.contains(&g.as_str()))
    {
        return None;
    }
    let mut via: Vec<String> = native.iter().filter(|g| !two_qubit(g)).cloned().collect();
    via.push("cx".to_string());
    Some(via)
}

/// The default pipeline (see [`compile`]), with permutation-aware block
/// resynthesis if `permute`.
fn pipeline(
    instructions: &[Instruction],
    native: &[String],
    permute: bool,
) -> Result<(Vec<Instruction>, Vec<usize>), String> {
    let mut lowering = Lowering::new(native.iter().cloned());
    let mut current = Vec::with_capacity(instructions.len());
    for inst in instructions {
        current.extend(lowering.lower(inst)?);
    }
    // Phase folding first, while rotations are still Z rotations between
    // CNOTs; each pass can expose work for the next.
    let cleanup = |current: Vec<Instruction>, lowering: &Lowering| {
        let current = crate::pauli_network::merge_rotations(current);
        let current = phase_folding(current);
        let current = merge_single_qubit_runs(current, &lowering.native, &lowering.bases_1q);
        commutative_cancellation(current)
    };
    let current = cleanup(current, &lowering);
    let (current, wires) = resynthesize_permuting(current, &mut lowering, permute);
    // Once more, without resynthesis: a second one never removed a gate on
    // the benchmarks, and resynthesis is most of the time.
    let current = cleanup(current, &lowering);
    let current = shift_through_symmetric_blocks(current, &lowering.bases_1q);
    let current = merge_single_qubit_runs(current, &lowering.native, &lowering.bases_1q);
    Ok((current, wires))
}

#[cfg(test)]
mod tests {
    use super::*;
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
            // Resynthesis may leave single-qubit gates for the next merge.
            let resynthesized = resynthesize_two_qubit_blocks(circuit.clone(), &mut lowering);
            let merge = |c: Vec<Instruction>| {
                merge_single_qubit_runs(c, &lowering.native, &lowering.bases_1q)
            };
            assert!(cost(&merge(resynthesized.clone())) <= cost(&merge(circuit.clone())));
            for out in [
                commutative_cancellation(circuit.clone()),
                resynthesized,
                merge(circuit.clone()),
            ] {
                assert!(dense_equal_up_to_phase(
                    &circuit_unitary(&out, 3),
                    &reference,
                    1e-8
                ));
            }
            assert!(cost(&commutative_cancellation(circuit.clone())) <= cost(&circuit));
            assert!(cost(&merge(circuit.clone())) <= cost(&circuit));
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
        let no_2q = vec![Instruction::gate("cx", vec![0, 1], vec![])];
        assert!(
            compile(no_2q, natives(&["rz", "sx", "swap"]))
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

    /// `circuit` followed by SWAPs that move every state back from `wires`.
    fn restored(mut circuit: Vec<Instruction>, mut wires: Vec<usize>) -> Vec<Instruction> {
        for q in 0..wires.len() {
            let w = wires[q];
            if w != q {
                circuit.push(Instruction::gate("swap", vec![w, q], vec![]));
                let p = wires.iter().position(|&x| x == q).expect("a permutation");
                wires[p] = w;
                wires[q] = q;
            }
        }
        circuit
    }

    #[test]
    fn a_cx_each_way_becomes_one_cx_and_a_relabelling() {
        let g = |name: &str, q: &[usize]| Instruction::gate(name, q.to_vec(), vec![]);
        let circuit = vec![
            g("h", &[0]),
            g("cx", &[0, 1]),
            g("cx", &[1, 0]),
            g("t", &[1]),
        ];
        let native = natives(&["rz", "sx", "x", "cx"]);
        let kept = compile(circuit.clone(), native.clone()).unwrap();
        let (out, wires) =
            compile_relabelling(circuit.clone(), native.clone(), Relabel::Always).unwrap();
        assert_eq!(cost(&kept).0, 2);
        assert_eq!(cost(&out).0, 1);
        assert_eq!(wires, vec![1, 0]);
        assert!(dense_equal_up_to_phase(
            &circuit_unitary(&restored(out, wires), 2),
            &circuit_unitary(&circuit, 2),
            1e-8
        ));
        // Unmeasured, the relabelling would show: IfMeasured keeps the CXs.
        let (unmeasured, wires) =
            compile_relabelling(circuit, native, Relabel::IfMeasured).unwrap();
        assert_eq!((cost(&unmeasured).0, wires), (2, vec![0, 1]));
    }

    /// Exponentials of random Pauli strings of weight 2 to 4 on `n` qubits,
    /// each as basis changes around a ladder of CXs, as Qiskit writes them.
    fn pauli_exponentials(rng: &mut Rng, n: usize, terms: usize) -> Vec<Instruction> {
        let g =
            |name: &str, q: &[usize], params: Vec<f64>| Instruction::gate(name, q.to_vec(), params);
        let mut circuit = Vec::new();
        for _ in 0..terms {
            let mut qubits: Vec<usize> = (0..n).collect();
            for i in 0..n {
                qubits.swap(i, i + rng.below(n - i));
            }
            qubits.truncate(2 + rng.below(3));
            let letters: Vec<usize> = qubits.iter().map(|_| rng.below(3)).collect();
            let into: Vec<Instruction> = qubits
                .iter()
                .zip(&letters)
                .flat_map(|(&q, &letter)| match letter {
                    0 => vec![g("h", &[q], vec![])],
                    1 => vec![g("sdg", &[q], vec![]), g("h", &[q], vec![])],
                    _ => vec![],
                })
                .collect();
            let ladder: Vec<Instruction> = qubits
                .windows(2)
                .map(|w| g("cx", &[w[0], w[1]], vec![]))
                .collect();
            let back = |inst: &Instruction| match inst.name.as_str() {
                "sdg" => g("s", &inst.qubits, vec![]),
                _ => inst.clone(),
            };
            circuit.extend(into.iter().cloned());
            circuit.extend(ladder.iter().cloned());
            let last = *qubits.last().expect("a qubit");
            circuit.push(g("rz", &[last], vec![rng.range(-3.0, 3.0)]));
            circuit.extend(ladder.iter().rev().cloned());
            circuit.extend(into.iter().rev().map(back));
        }
        circuit
    }

    #[test]
    fn pauli_networks_are_kept_only_when_they_cost_no_more() {
        let mut rng = Rng::new(23);
        let native = natives(&["rz", "sx", "x", "cx"]);
        let mut kept = 0;
        for _ in 0..12 {
            let n = 5;
            let circuit = pauli_exponentials(&mut rng, n, 6);
            let (out, wires) =
                compile_relabelling(circuit.clone(), native.clone(), Relabel::Always).unwrap();
            assert!(dense_equal_up_to_phase(
                &circuit_unitary(&restored(out.clone(), wires.clone()), n),
                &circuit_unitary(&circuit, n),
                1e-7
            ));
            // No more two-qubit gates, nor gates in all, than without it.
            let (direct, _) = pipeline(&circuit, &native, true).unwrap();
            let ((two, all), (direct_two, direct_all)) = (cost(&out), cost(&direct));
            assert!(two <= direct_two && all <= direct_all);
            kept += usize::from(two < direct_two);
            // The seeded trials run side by side, but the result is always the same.
            let again = compile_relabelling(circuit, native.clone(), Relabel::Always).unwrap();
            assert_eq!((out, wires), again);
        }
        assert!(kept >= 6, "the network was kept {kept} times in 12");
    }

    #[test]
    fn relabelled_compilations_agree_up_to_the_permutation() {
        let mut rng = Rng::new(17);
        let native = natives(&["rz", "sx", "x", "cx"]);
        for _ in 0..60 {
            let circuit = random_circuit(&mut rng, 3, 25);
            let (out, wires) =
                compile_relabelling(circuit.clone(), native.clone(), Relabel::Always).unwrap();
            assert!(dense_equal_up_to_phase(
                &circuit_unitary(&restored(out, wires), 3),
                &circuit_unitary(&circuit, 3),
                1e-7
            ));
        }
    }

    #[test]
    fn merged_cost_counts_what_the_merge_leaves() {
        let mut rng = Rng::new(9);
        for native in [
            &["rz", "sx", "x", "cx"][..],
            &["u3", "cx"],
            &["rx", "ry", "rxx"],
        ] {
            let lowering = Lowering::new(natives(native));
            for _ in 0..50 {
                let circuit = random_circuit(&mut rng, 3, 30);
                let merged =
                    merge_single_qubit_runs(circuit.clone(), &lowering.native, &lowering.bases_1q);
                assert_eq!(
                    merged_cost(&circuit, &lowering.native, &lowering.bases_1q),
                    cost(&merged)
                );
            }
        }
    }
}
