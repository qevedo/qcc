//! Gathering two-qubit blocks through commuting gates.
//!
//! Block resynthesis (`passes.rs`) takes maximal runs of gates on one pair of
//! qubits as they stand in the list. Two gates on the same pair that a gate
//! on one of their qubits and a third qubit separates fall in different
//! blocks, even when that gate commutes with everything between it and the
//! end of the block: `cx a,b; cx a,c; cx b,a` is two blocks, but `cx a,c`
//! commutes with `cx a,b`, so it is also `cx a,c; cx a,b; cx b,a`, one block
//! of two CXs that a SWAP away needs only one. Clifford simplification in
//! tket (Sivarajah et al., arXiv:2003.10611; Fagan and Duncan,
//! arXiv:1901.10114) finds such pairs by moving each CX's Pauli
//! interactions through commuting gates; this pass does the same for any
//! two-qubit gate, and lets the resynthesis do the rest.
//!
//! For each two-qubit gate `g` on a pair `{a, b}`, the pass looks back for the
//! previous two-qubit gate `h` on the same pair. The two-qubit gates between
//! them that touch `a` or `b` (the blockers, at most [`MAX_BLOCKERS`]) are
//! moved out: some to before `h`, the others to after `g`, or further, into
//! the runs of gates on `{a, b}` that end at `h` and start at `g`, wherever
//! the cut leaves the fewest CX-like gates in the blocks it makes. A blocker
//! moves only if it commutes, exactly, with the product of the pair's gates
//! it crosses (one 8×8 check: single-qubit gates on the shared qubit count
//! together, not one by one). The gates on its third qubit that it depends
//! on (or that depend on it) move with it; they share no qubit with the
//! pair, so for them only the order of the list changes, as in the circuit's
//! DAG nothing does. The unitary is unchanged.
//!
//! A move also breaks the blocks the blockers were in, so the pass only
//! makes it when the gathered blocks need fewer CX-like gates than the two
//! runs they come from (with `permute`, counting a block that a SWAP away
//! needs fewer as the cheaper of the two, as the resynthesis will elide that
//! SWAP). Whether the whole circuit gets cheaper is for the pipeline to
//! judge (`pipeline` in `passes.rs`).
//!
//! Where the blockers do not commute with the blocks, part of a block may
//! still commute with them once moved through the single-qubit gates
//! between: [`carry_through`] carries the end of a block forward, or the
//! start of the next one back, so that the two blocks join.

use crate::gates::{is_known_gate, matrix_1q, matrix_2q};
use crate::ir::Instruction;
use crate::linalg::{C, M2, M4, kron};
use crate::one_qubit::{Basis1q, synthesize_1q};
use crate::passes::pair_matrix;
use crate::two_qubit::cx_count;

/// How many blockers a pair of gates may have between them: more are rarely
/// all movable, and each split of them into left and right is tried.
const MAX_BLOCKERS: usize = 3;

/// How many places in each run a gathered block may start or end at, the
/// nearest to the blockers first.
const CUTS: usize = 8;

/// How many instructions on the pair's qubits to look through, each way.
const WINDOW: usize = 64;

/// A 3-qubit operator; qubit 0 is the most significant.
type M8 = [[C; 8]; 8];

/// `instructions` with the gates that keep two-qubit blocks on the same pair
/// apart moved out of the way, where they commute and the joined block
/// needs fewer two-qubit gates (see the module documentation).
pub fn gather_blocks(instructions: Vec<Instruction>, permute: bool) -> Vec<Instruction> {
    let mut insts = instructions;
    // Each gathering is judged on its own, so in principle a blocker moved
    // out of one block could later gather its own block back across it:
    // a budget of one gathering per instruction rules out cycling.
    let mut budget = insts.len();
    let mut i = 0;
    while i < insts.len() {
        i = match gather_at(&mut insts, i, permute) {
            Some(next) if budget > 0 => {
                budget -= 1;
                next
            }
            Some(_) => break,
            None => i + 1,
        };
    }
    insts
}

fn two_qubit_gate(inst: &Instruction) -> bool {
    inst.qubits.len() == 2 && inst.is_gate() && is_known_gate(&inst.name)
}

fn touches(inst: &Instruction, pair: [usize; 2]) -> bool {
    inst.qubits.iter().any(|q| pair.contains(q))
}

/// A gate acting on the pair's qubits only.
fn on_pair(inst: &Instruction, pair: [usize; 2]) -> bool {
    inst.is_gate()
        && is_known_gate(&inst.name)
        && !inst.qubits.is_empty()
        && inst.qubits.iter().all(|q| pair.contains(q))
}

/// Gather the block of the two-qubit gate at `i` with the previous one on its
/// pair, if that pays; returns where to continue.
fn gather_at(insts: &mut Vec<Instruction>, i: usize, permute: bool) -> Option<usize> {
    if !two_qubit_gate(&insts[i]) {
        return None;
    }
    let pair = [insts[i].qubits[0], insts[i].qubits[1]];
    // Back to the previous two-qubit gate on the pair, through blockers only.
    let mut blockers = Vec::new();
    let mut previous = None;
    let mut seen = 0;
    for k in (0..i).rev() {
        let x = &insts[k];
        if !touches(x, pair) {
            continue;
        }
        seen += 1;
        if seen > WINDOW {
            return None;
        }
        if on_pair(x, pair) {
            if x.qubits.len() == 2 {
                previous = Some(k);
                break;
            }
        } else if two_qubit_gate(x) {
            blockers.push(k);
            if blockers.len() > MAX_BLOCKERS {
                return None;
            }
        } else {
            return None;
        }
    }
    let h = previous?;
    if blockers.is_empty() {
        return None;
    }
    blockers.reverse();
    // The runs on the pair that end at h and start at i.
    let run_end = |range: &mut dyn Iterator<Item = usize>, insts: &[Instruction]| {
        let mut end = None;
        let mut seen = 0;
        for k in range {
            let x = &insts[k];
            if !touches(x, pair) {
                continue;
            }
            seen += 1;
            if !on_pair(x, pair) || seen > WINDOW {
                break;
            }
            end = Some(k);
        }
        end
    };
    let start = run_end(&mut (0..h).rev(), insts).unwrap_or(h);
    let end = run_end(&mut (i + 1..insts.len()), insts).unwrap_or(i);
    let is_member = |k: usize| touches(&insts[k], pair) && !blockers.contains(&k);
    // Each cut's cost is that of three ranges, and the cuts share their
    // outer ranges: each range's cost is worked out once (its unitary and
    // CX count are the same every time, so remembering them changes no
    // choice).
    let costs = std::cell::RefCell::new(std::collections::HashMap::new());
    let cost = |range: std::ops::Range<usize>| -> Option<usize> {
        let key = (range.start, range.end);
        if let Some(&known) = costs.borrow().get(&key) {
            return known;
        }
        let members = range.filter(|&k| is_member(k)).map(|k| &insts[k]);
        let known = pair_matrix(pair, members).map(|m| cx_cost(&m, permute));
        costs.borrow_mut().insert(key, known);
        known
    };
    // The blockers may stop anywhere in the two runs: the gathered block
    // starts before member `first` of the run that ends at h and ends after
    // member `last` of the run that starts at i, and the rest of each run
    // stays a block of its own. Blockers moved only across h and i have
    // the fewest gates to commute with; across the whole runs, the joined
    // block is largest. The region with the first `split` blockers moved
    // left and the others right, for the first split that can move. The
    // region is the circuit's indices from start to end, reordered: moving
    // gates changes only their order, so no instruction is copied until the
    // chosen order replaces the region.
    let region: Vec<usize> = (start..end + 1).collect();
    let moved = |first: usize, last: usize| -> Option<Vec<usize>> {
        (0..=blockers.len()).rev().find_map(|split| {
            let mut work = region.clone();
            for &b in &blockers[..split] {
                if !move_blocker(&mut work, insts, b, first, pair, true) {
                    return None;
                }
            }
            for &b in blockers[split..].iter().rev() {
                if !move_blocker(&mut work, insts, b, last, pair, false) {
                    return None;
                }
            }
            Some(work)
        })
    };
    // Blockers that cannot pass even h and i alone rarely pass more (on
    // Quantum Volume's generic blocks, never): stop before any synthesis or
    // copying. The gates they carry do not change what they cross.
    let crossing = |p: usize, range: std::ops::Range<usize>| {
        let crossed = range.filter(|&k| is_member(k)).map(|k| &insts[k]).collect();
        crosses(&insts[p], crossed, pair)
    };
    (0..=blockers.len()).find(|&split| {
        blockers[..split].iter().all(|&p| crossing(p, h..p))
            && blockers[split..].iter().all(|&p| crossing(p, p + 1..i + 1))
    })?;
    let current = cost(start..blockers[0])? + cost(blockers[blockers.len() - 1] + 1..end + 1)?;
    let firsts: Vec<usize> = (start..h + 1).filter(|&k| is_member(k)).collect();
    let lasts: Vec<usize> = (i..end + 1).filter(|&k| is_member(k)).collect();
    // The cheapest cut, the widest on ties.
    let mut best: Option<(usize, Vec<usize>)> = None;
    for &first in firsts.iter().skip(firsts.len().saturating_sub(CUTS)) {
        for &last in lasts.iter().take(CUTS).rev() {
            let total = cost(start..first)? + cost(first..last + 1)? + cost(last + 1..end + 1)?;
            if total >= best.as_ref().map_or(current, |b| b.0) {
                continue;
            }
            if let Some(work) = moved(first, last) {
                best = Some((total, work));
            }
        }
    }
    let (_, work) = best?;
    let g = start + work.iter().position(|&k| k == i).expect("the gate itself");
    let gathered: Vec<Instruction> = work.iter().map(|&k| insts[k].clone()).collect();
    insts.splice(start..end + 1, gathered);
    Some(g + 1)
}

/// Move the blocker with index `blocker` in `work` (indices into `insts`,
/// in their current order) to just before the instruction `anchor`, or
/// just after it if not `leftward`. It must commute with the product of the
/// gates on `pair` it crosses. The gates on its other qubit that it depends
/// on (or that depend on it, moving right), and in turn on theirs, move with
/// it: they share no qubit with the gates they cross, so only the order of
/// the list changes for them, not the circuit.
fn move_blocker(
    work: &mut Vec<usize>,
    insts: &[Instruction],
    blocker: usize,
    anchor: usize,
    pair: [usize; 2],
    leftward: bool,
) -> bool {
    let position = |work: &[usize], tag: usize| {
        work.iter()
            .position(|&k| k == tag)
            .expect("an index of the region")
    };
    let (p, a) = (position(work, blocker), position(work, anchor));
    let gate = &insts[blocker];
    let mut reach: Vec<usize> = gate
        .qubits
        .iter()
        .filter(|q| !pair.contains(q))
        .copied()
        .collect();
    let mut carried = vec![p];
    let mut crossed = Vec::new();
    let range: Vec<usize> = if leftward {
        (a..p).rev().collect()
    } else {
        (p + 1..a + 1).collect()
    };
    for k in range {
        let x = &insts[work[k]];
        if touches(x, pair) {
            if !on_pair(x, pair) {
                return false;
            }
            crossed.push(x);
        } else if x.qubits.iter().any(|q| reach.contains(q)) {
            if !x.is_gate() {
                return false;
            }
            reach.extend(x.qubits.iter().copied());
            carried.push(k);
        }
    }
    if leftward {
        crossed.reverse();
    }
    if !crosses(gate, crossed, pair) {
        return false;
    }
    carried.sort_unstable();
    let group: Vec<usize> = carried.iter().map(|&k| work[k]).collect();
    let mut rest: Vec<usize> = work
        .iter()
        .enumerate()
        .filter(|(k, _)| !carried.contains(k))
        .map(|(_, &x)| x)
        .collect();
    let at = position(&rest, anchor) + usize::from(!leftward);
    rest.splice(at..at, group);
    *work = rest;
    true
}

/// How many CX-like gates the resynthesis will spend on a block: with
/// `permute`, it may synthesize SWAP·U instead and elide the SWAP.
fn cx_cost(m: &M4, permute: bool) -> usize {
    let direct = cx_count(m);
    if permute {
        let swap = matrix_2q("swap", &[]).expect("a standard gate");
        direct.min(cx_count(&(swap * *m)))
    } else {
        direct
    }
}

/// Whether `blocker`, which touches one qubit of `pair`, commutes with the
/// product of `crossed`, gates on `pair` in time order.
fn crosses(blocker: &Instruction, crossed: Vec<&Instruction>, pair: [usize; 2]) -> bool {
    let Some(shared) = blocker.qubits.iter().find(|q| pair.contains(q)).copied() else {
        return false;
    };
    // Only the gates acting on the shared qubit can fail to commute.
    if !crossed.iter().any(|x| x.qubits.contains(&shared)) {
        return true;
    }
    let Some(s) = pair_matrix(pair, crossed) else {
        return false;
    };
    commutes_on_pair(blocker, &s, pair)
}

/// Whether `blocker`, a two-qubit gate on one qubit of `pair` and another
/// qubit, commutes with `s`, an operator on `pair` (`pair[0]` the most
/// significant qubit).
fn commutes_on_pair(blocker: &Instruction, s: &M4, pair: [usize; 2]) -> bool {
    let Some(shared) = blocker.qubits.iter().find(|q| pair.contains(q)).copied() else {
        return false;
    };
    let Some(other) = blocker.qubits.iter().find(|&&q| q != shared).copied() else {
        return false;
    };
    let unshared = if pair[0] == shared { pair[1] } else { pair[0] };
    let qubits = [shared, unshared, other];
    let s = embed(s, [position(&qubits, pair[0]), position(&qubits, pair[1])]);
    on3(blocker, &qubits).is_some_and(|b| commutes(&b, &s))
}

fn position(qubits: &[usize], q: usize) -> usize {
    qubits
        .iter()
        .position(|&x| x == q)
        .expect("a qubit of the three")
}

/// A one- or two-qubit gate as an operator on `qubits` (three of them).
fn on3(inst: &Instruction, qubits: &[usize]) -> Option<M8> {
    match *inst.qubits.as_slice() {
        [q] => {
            let g = matrix_1q(&inst.name, &inst.params)?;
            let mut m = M4::IDENTITY;
            // As a two-qubit operator on (q, a qubit it does not touch).
            for r in 0..4 {
                for c in 0..4 {
                    m.0[r][c] = if r % 2 == c % 2 {
                        g.0[r / 2][c / 2]
                    } else {
                        C::new(0.0, 0.0)
                    };
                }
            }
            let p = position(qubits, q);
            Some(embed(&m, [p, (p + 1) % 3]))
        }
        [q0, q1] => {
            let g = matrix_2q(&inst.name, &inst.params)?;
            Some(embed(&g, [position(qubits, q0), position(qubits, q1)]))
        }
        _ => None,
    }
}

/// A two-qubit operator acting on positions `at` (most significant first)
/// of three qubits.
fn embed(m: &M4, at: [usize; 2]) -> M8 {
    let bit = |index: usize, p: usize| (index >> (2 - p)) & 1;
    let rest = 3 - at[0] - at[1];
    let mut out = [[C::new(0.0, 0.0); 8]; 8];
    for i in 0..8 {
        for j in 0..8 {
            if bit(i, rest) == bit(j, rest) {
                let r = (bit(i, at[0]) << 1) | bit(i, at[1]);
                let c = (bit(j, at[0]) << 1) | bit(j, at[1]);
                out[i][j] = m.0[r][c];
            }
        }
    }
    out
}

fn mul8(a: &M8, b: &M8) -> M8 {
    let mut out = [[C::new(0.0, 0.0); 8]; 8];
    for i in 0..8 {
        for k in 0..8 {
            if a[i][k] == C::new(0.0, 0.0) {
                continue;
            }
            for j in 0..8 {
                out[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    out
}

fn commutes(a: &M8, b: &M8) -> bool {
    let (ab, ba) = (mul8(a, b), mul8(b, a));
    (0..8).all(|i| (0..8).all(|j| (ab[i][j] - ba[i][j]).norm() < 1e-9))
}

// -- carrying gates of a block to the next block on its pair ------------------------------

/// How many two-qubit gates at the end (or start) of a block may be carried
/// to the next (or previous) block on its pair at once.
const CARRIED: usize = 3;

/// `instructions` with the last two-qubit gates of a block carried forward
/// into the next block on the same pair, or the first ones of a block
/// carried back into the previous one, where that leaves fewer CX-like
/// gates in the two blocks.
///
/// Gathering ([`gather_blocks`]) moves the gates between two blocks on a
/// pair out of the way, whole: it cannot when they do not commute with the
/// blocks. But a part of a block may commute with them, once moved through
/// the single-qubit gates between. In `cx a,b; h b; cx b,c; cx a,b` the
/// blocker `cx b,c` commutes with neither `cx a,b` nor `h b`, yet
/// `h b · cx a,b = cz a,b · h b`, and the CZ, diagonal, commutes with the
/// control of `cx b,c`: the circuit is `h b; cx b,c; (cx a,b · cz a,b)`,
/// two CXs instead of three. ZX-calculus Clifford simplification (Fagan
/// and Duncan, arXiv:1901.10114) makes rewrites of this kind on Clifford
/// circuits; here any gates of a block may move, with exact checks:
///
/// - A carried part `D`, the gates of a block from one of its last
///   [`CARRIED`] two-qubit gates to its end, crosses each single-qubit gate
///   `g` after it on the pair by becoming `g D g†`, and each blocker only if
///   it commutes with `D` as it stands there (an 8×8 check).
/// - At the next block it is put down as the gates it was, between the
///   inverse of the single-qubit gates it crossed and those gates again
///   (`R† D R` in time order), for the resynthesis to merge.
/// - The same backwards: the first gates of a block, carried back to just
///   after the previous one.
///
/// A move is made only when the two blocks it changes need fewer CX-like
/// gates than before (with `permute`, as [`gather_blocks`] counts them).
/// The single-qubit gates the merged blocks need can grow, so whether the
/// circuit gets cheaper is for the pipeline to judge.
pub fn carry_through(
    instructions: Vec<Instruction>,
    permute: bool,
    bases: &[Basis1q],
) -> Vec<Instruction> {
    if bases.is_empty() {
        return instructions;
    }
    let mut insts = instructions;
    // Each move lowers the CX count of the blocks it changes, so the pass
    // ends; the budget is a guard all the same.
    let mut budget = insts.len();
    let mut i = 0;
    while i < insts.len() {
        i = match carry_at(&mut insts, i, permute, bases) {
            Some(next) if budget > 0 => {
                budget -= 1;
                next
            }
            Some(_) => break,
            None => i + 1,
        };
    }
    insts
}

/// A single-qubit gate on one qubit of `pair` as an operator on the pair.
fn on_pair_1q(inst: &Instruction, pair: [usize; 2]) -> Option<M4> {
    let g = matrix_1q(&inst.name, &inst.params)?;
    Some(if inst.qubits[0] == pair[0] {
        kron(&g, &M2::IDENTITY)
    } else {
        kron(&M2::IDENTITY, &g)
    })
}

/// Carry gates between the block of the two-qubit gate at `i` and the next
/// block on its pair, if that pays; returns where to continue.
fn carry_at(
    insts: &mut Vec<Instruction>,
    i: usize,
    permute: bool,
    bases: &[Basis1q],
) -> Option<usize> {
    if !two_qubit_gate(&insts[i]) {
        return None;
    }
    let pair = [insts[i].qubits[0], insts[i].qubits[1]];
    // The gap: what touches the pair between `i` and the next two-qubit
    // gate on it, `j`. Single-qubit gates and blockers only.
    let mut gap = Vec::new();
    let mut blockers = 0;
    let mut j = None;
    for k in i + 1..insts.len() {
        let x = &insts[k];
        if !touches(x, pair) {
            continue;
        }
        if gap.len() >= WINDOW {
            return None;
        }
        if on_pair(x, pair) {
            if x.qubits.len() == 2 {
                j = Some(k);
                break;
            }
            gap.push(k);
        } else if two_qubit_gate(x) {
            blockers += 1;
            if blockers > MAX_BLOCKERS {
                return None;
            }
            gap.push(k);
        } else {
            return None;
        }
    }
    let j = j?;
    if blockers == 0 {
        return None;
    }
    // The blocks: the run on the pair that ends at i, and the one that
    // starts at j.
    let run = |range: &mut dyn Iterator<Item = usize>| {
        let mut members = Vec::new();
        for k in range {
            let x = &insts[k];
            if !touches(x, pair) {
                continue;
            }
            if !on_pair(x, pair) || members.len() >= WINDOW {
                break;
            }
            members.push(k);
        }
        members
    };
    let mut first: Vec<usize> = vec![i];
    first.extend(run(&mut (0..i).rev()));
    first.reverse();
    let mut second: Vec<usize> = vec![j];
    second.extend(run(&mut (j + 1..insts.len())));
    let matrix = |members: &[usize]| pair_matrix(pair, members.iter().map(|&k| &insts[k]));
    // The gap as steps: each single-qubit gate as an operator on the pair,
    // each blocker as None. And the single-qubit gates on each qubit of the
    // pair, multiplied in time order.
    let mut steps: Vec<(usize, Option<M4>)> = Vec::with_capacity(gap.len());
    let mut crossed = [M2::IDENTITY; 2];
    for &k in &gap {
        let x = &insts[k];
        if x.qubits.len() == 1 {
            let g = matrix_1q(&x.name, &x.params)?;
            let side = usize::from(x.qubits[0] == pair[1]);
            crossed[side] = g * crossed[side];
            steps.push((k, Some(on_pair_1q(x, pair)?)));
        } else {
            steps.push((k, None));
        }
    }
    // `d` moved through the gap, forward or backward, as it ends up; None if
    // a blocker on the way does not commute with it as it stands there.
    let carry = |mut d: M4, forward: bool| -> Option<M4> {
        let mut cross = |&(k, e): &(usize, Option<M4>)| -> Option<()> {
            match e {
                Some(e) if forward => d = e * d * e.adjoint(),
                Some(e) => d = e.adjoint() * d * e,
                None if commutes_on_pair(&insts[k], &d, pair) => {}
                None => return None,
            }
            Some(())
        };
        if forward {
            steps.iter().try_for_each(&mut cross)?;
        } else {
            steps.iter().rev().try_for_each(&mut cross)?;
        }
        Some(d)
    };
    let two_qubit = |members: &[usize]| -> Vec<usize> {
        (0..members.len())
            .filter(|&p| insts[members[p]].qubits.len() == 2)
            .collect()
    };
    // The parts that can be carried, as (forward, members moved, members
    // left in their block, the moved part as it ends up). Forward: the end
    // of the first block, from one of its last two-qubit gates, joins the
    // second. Backward: the start of the second block joins the first.
    let mut movable: Vec<(bool, &[usize], &[usize], M4)> = Vec::new();
    for &cut in two_qubit(&first).iter().rev().take(CARRIED) {
        let (left, moved) = first.split_at(cut);
        if let Some(d) = matrix(moved).and_then(|d| carry(d, true)) {
            movable.push((true, moved, left, d));
        }
    }
    for &cut in two_qubit(&second).iter().take(CARRIED) {
        let (moved, right) = second.split_at(cut + 1);
        if let Some(d) = matrix(moved).and_then(|d| carry(d, false)) {
            movable.push((false, moved, right, d));
        }
    }
    // Most gates cannot be carried at all (on Quantum Volume's generic
    // blocks, almost none): stop before any KAK decomposition.
    if movable.is_empty() {
        return None;
    }
    let (m_first, m_second) = (matrix(&first)?, matrix(&second)?);
    let current = cx_cost(&m_first, permute) + cx_cost(&m_second, permute);
    // The fewest CX-like gates, then the fewest gates moved.
    let mut best: Option<(usize, bool, Vec<usize>)> = None;
    for (forward, moved, rest, d) in movable {
        let Some(rest) = matrix(rest) else {
            continue;
        };
        let joined = if forward { m_second * d } else { d * m_first };
        let total = cx_cost(&rest, permute) + cx_cost(&joined, permute);
        let better = match &best {
            None => total < current,
            Some((t, _, m)) => (total, moved.len()) < (*t, m.len()),
        };
        if better {
            best = Some((total, forward, moved.to_vec()));
        }
    }
    let (_, forward, moved) = best?;
    // The moved gates, between the gap's single-qubit gates undone and
    // redone: forward, R† D R in time order; backward, R D R†.
    let emit_all = |ms: [M2; 2]| -> Vec<Instruction> {
        pair.iter()
            .zip(ms)
            .flat_map(|(&q, m)| emit_1q(q, &m, bases))
            .collect()
    };
    let undone = [crossed[0].adjoint(), crossed[1].adjoint()];
    let (before, after) = if forward {
        (emit_all(undone), emit_all(crossed))
    } else {
        (emit_all(crossed), emit_all(undone))
    };
    let mut put: Vec<Instruction> = before;
    put.extend(moved.iter().map(|&k| insts[k].clone()));
    put.extend(after);
    // Rebuild the list: the moved gates leave, and come back just before j
    // (forward) or just after i (backward).
    let at = if forward { j } else { i + 1 };
    let mut out = Vec::with_capacity(insts.len() + put.len());
    for (k, inst) in insts.drain(..).enumerate() {
        if k == at {
            out.append(&mut put);
        }
        if !moved.contains(&k) {
            out.push(inst);
        }
    }
    *insts = out;
    Some(first[0])
}

fn emit_1q(qubit: usize, m: &M2, bases: &[Basis1q]) -> Vec<Instruction> {
    synthesize_1q(m, bases)
        .into_iter()
        .map(|op| Instruction::gate(&op.name, [qubit], op.params))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testing::{Rng, circuit_unitary, dense_equal_up_to_phase};

    fn gate(name: &str, qubits: &[usize], params: &[f64]) -> Instruction {
        Instruction::gate(name, qubits, params)
    }

    #[test]
    fn a_commuting_cx_moves_out_of_a_block() {
        // cx 0,2 commutes with cx 0,1, and the two CXs on {0, 1} are a SWAP
        // and one CX: with permute, they are worth gathering.
        let circuit = vec![
            gate("cx", &[0, 1], &[]),
            gate("cx", &[0, 2], &[]),
            gate("cx", &[1, 0], &[]),
        ];
        let out = gather_blocks(circuit.clone(), true);
        // cx 0,2 does not commute with cx 1,0 (its target is the control),
        // so it moves left.
        assert_eq!(out[0].qubits, vec![0, 2]);
        assert!(dense_equal_up_to_phase(
            &circuit_unitary(&circuit, 3),
            &circuit_unitary(&out, 3),
            1e-9
        ));
        // Without permute, two CXs are no worse than one CX each side.
        assert_eq!(gather_blocks(circuit.clone(), false), circuit);
    }

    #[test]
    fn a_blocker_that_does_not_commute_stays() {
        // cx 2,0 targets qubit 0, which cx 0,1 controls: they do not commute
        // either way.
        let circuit = vec![
            gate("cx", &[0, 1], &[]),
            gate("h", &[1], &[]),
            gate("cx", &[2, 0], &[]),
            gate("cx", &[0, 1], &[]),
            gate("h", &[1], &[]),
            gate("cx", &[0, 1], &[]),
        ];
        assert_eq!(gather_blocks(circuit.clone(), true), circuit);
    }

    #[test]
    fn single_qubit_gates_count_together() {
        // cx 0,2 cannot move left past h 0, nor right past either x on its
        // own, but the product of the gates on {0, 1} it passes, x x cx 0,1,
        // commutes with it. The joined block is h 0 alone.
        let circuit = vec![
            gate("h", &[0], &[]),
            gate("cx", &[0, 1], &[]),
            gate("cx", &[0, 2], &[]),
            gate("x", &[0], &[]),
            gate("x", &[0], &[]),
            gate("cx", &[0, 1], &[]),
        ];
        let out = gather_blocks(circuit.clone(), false);
        assert_eq!(out.last().unwrap().qubits, vec![0, 2]);
        assert!(dense_equal_up_to_phase(
            &circuit_unitary(&circuit, 3),
            &circuit_unitary(&out, 3),
            1e-9
        ));
    }

    #[test]
    fn gathering_keeps_random_circuits_equivalent() {
        // Clifford+T-like circuits on few qubits, where gathering happens
        // often; every output must have the input's unitary.
        let mut rng = Rng::new(17);
        let names = ["h", "s", "t", "sx", "x", "rz"];
        let mut gathered = 0;
        for trial in 0..300 {
            let n = 3 + trial % 3;
            let mut circuit = Vec::new();
            for _ in 0..30 {
                if rng.uniform() < 0.5 {
                    let a = rng.below(n);
                    let b = (a + 1 + rng.below(n - 1)) % n;
                    circuit.push(gate(["cx", "cz"][rng.below(2)], &[a, b], &[]));
                } else {
                    let name = names[rng.below(names.len())];
                    let params = if name == "rz" {
                        vec![rng.uniform()]
                    } else {
                        vec![]
                    };
                    circuit.push(gate(name, &[rng.below(n)], &params));
                }
            }
            for permute in [false, true] {
                let out = gather_blocks(circuit.clone(), permute);
                gathered += usize::from(out != circuit);
                assert!(dense_equal_up_to_phase(
                    &circuit_unitary(&circuit, n),
                    &circuit_unitary(&out, n),
                    1e-9
                ));
            }
        }
        assert!(gathered > 20, "only {gathered} circuits changed");
    }

    fn ibm_bases() -> Vec<Basis1q> {
        let native = ["rz", "sx", "x"].iter().map(|g| g.to_string()).collect();
        crate::one_qubit::bases_for(&native)
    }

    /// The CX-like gates the blocks on `pair` need, in list order.
    fn pair_costs(circuit: &[Instruction], pair: [usize; 2]) -> Vec<usize> {
        let mut costs = Vec::new();
        let mut run: Vec<&Instruction> = Vec::new();
        for inst in circuit {
            if !touches(inst, pair) {
                continue;
            }
            if on_pair(inst, pair) {
                run.push(inst);
            } else {
                if run.iter().any(|i| i.qubits.len() == 2) {
                    costs.push(cx_count(&pair_matrix(pair, run.iter().copied()).unwrap()));
                }
                run.clear();
            }
        }
        if run.iter().any(|i| i.qubits.len() == 2) {
            costs.push(cx_count(&pair_matrix(pair, run.iter().copied()).unwrap()));
        }
        costs
    }

    #[test]
    fn a_cx_carried_through_a_hadamard_joins_the_next_block() {
        // h 1 · cx 0,1 is cz 0,1 · h 1, and the CZ commutes with the control
        // of cx 1,2: the first CX joins the last, and cx·cz on {0, 1} is one
        // CX. Neither gate commutes with cx 1,2 as it stands, so gathering
        // cannot move it.
        let circuit = vec![
            gate("cx", &[0, 1], &[]),
            gate("h", &[1], &[]),
            gate("cx", &[1, 2], &[]),
            gate("cx", &[0, 1], &[]),
        ];
        assert_eq!(gather_blocks(circuit.clone(), false), circuit);
        let out = carry_through(circuit.clone(), false, &ibm_bases());
        assert!(dense_equal_up_to_phase(
            &circuit_unitary(&circuit, 3),
            &circuit_unitary(&out, 3),
            1e-9
        ));
        assert_eq!(pair_costs(&circuit, [0, 1]), vec![1, 1]);
        assert_eq!(pair_costs(&out, [0, 1]), vec![1]);
        // The blocker now comes before every gate on {0, 1}.
        let blocker = out.iter().position(|i| i.qubits == [1, 2]).unwrap();
        assert!(out[..blocker].iter().all(|i| i.qubits.len() == 1));
    }

    #[test]
    fn a_cx_carried_back_through_a_hadamard_joins_the_previous_block() {
        // The mirror image: cx 0,1 · h 1 is h 1 · cz 0,1.
        let circuit = vec![
            gate("cx", &[0, 1], &[]),
            gate("cx", &[1, 2], &[]),
            gate("h", &[1], &[]),
            gate("cx", &[0, 1], &[]),
        ];
        let out = carry_through(circuit.clone(), false, &ibm_bases());
        assert!(dense_equal_up_to_phase(
            &circuit_unitary(&circuit, 3),
            &circuit_unitary(&out, 3),
            1e-9
        ));
        assert_eq!(pair_costs(&out, [0, 1]), vec![1]);
        let blocker = out.iter().position(|i| i.qubits == [1, 2]).unwrap();
        assert!(out[blocker + 1..].iter().all(|i| i.qubits.len() == 1));
    }

    #[test]
    fn nothing_is_carried_through_a_gate_it_does_not_commute_with() {
        // cx 0,1 targets qubit 1, which cx 1,2 controls, and no single-qubit
        // gate stands between to change that.
        let circuit = vec![
            gate("cx", &[0, 1], &[]),
            gate("cx", &[1, 2], &[]),
            gate("cx", &[0, 1], &[]),
        ];
        assert_eq!(carry_through(circuit.clone(), false, &ibm_bases()), circuit);
        // Wherever gates move, the blocks on the pair need fewer CXs.
        let circuit = vec![
            gate("cx", &[0, 1], &[]),
            gate("h", &[1], &[]),
            gate("cx", &[1, 2], &[]),
            gate("h", &[0], &[]),
            gate("cz", &[1, 0], &[]),
            gate("rz", &[1], &[0.3]),
            gate("cx", &[1, 0], &[]),
        ];
        let out = carry_through(circuit.clone(), false, &ibm_bases());
        let before: usize = pair_costs(&circuit, [0, 1]).iter().sum();
        let after: usize = pair_costs(&out, [0, 1]).iter().sum();
        assert!(out == circuit || after < before);
    }

    #[test]
    fn carrying_keeps_random_circuits_equivalent() {
        // As for gathering: Clifford+T-like circuits on few qubits, where
        // gates can often be carried; every output must have the input's
        // unitary.
        let mut rng = Rng::new(23);
        let names = ["h", "s", "t", "sx", "x", "rz"];
        let bases = ibm_bases();
        let mut carried = 0;
        for trial in 0..300 {
            let n = 3 + trial % 3;
            let mut circuit = Vec::new();
            for _ in 0..30 {
                if rng.uniform() < 0.5 {
                    let a = rng.below(n);
                    let b = (a + 1 + rng.below(n - 1)) % n;
                    circuit.push(gate(["cx", "cz"][rng.below(2)], &[a, b], &[]));
                } else {
                    let name = names[rng.below(names.len())];
                    let params = if name == "rz" {
                        vec![rng.uniform()]
                    } else {
                        vec![]
                    };
                    circuit.push(gate(name, &[rng.below(n)], &params));
                }
            }
            for permute in [false, true] {
                let out = carry_through(circuit.clone(), permute, &bases);
                carried += usize::from(out != circuit);
                assert!(dense_equal_up_to_phase(
                    &circuit_unitary(&circuit, n),
                    &circuit_unitary(&out, n),
                    1e-9
                ));
            }
        }
        assert!(carried > 20, "only {carried} circuits changed");
    }
}
