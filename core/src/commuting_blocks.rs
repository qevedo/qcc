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

use crate::gates::{is_known_gate, matrix_1q, matrix_2q};
use crate::ir::Instruction;
use crate::linalg::{C, M4};
use crate::passes::Block;
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
    let members = |range: std::ops::Range<usize>| -> Vec<Instruction> {
        range
            .filter(|&k| is_member(k))
            .map(|k| insts[k].clone())
            .collect()
    };
    let cost = |range: std::ops::Range<usize>| -> Option<usize> {
        let instructions = members(range);
        let m = Block { pair, instructions }.matrix()?;
        Some(cx_cost(&m, permute))
    };
    // The blockers may stop anywhere in the two runs: the gathered block
    // starts before member `first` of the run that ends at h and ends after
    // member `last` of the run that starts at i, and the rest of each run
    // stays a block of its own. Blockers moved only across h and i have
    // the fewest gates to commute with; across the whole runs, the joined
    // block is largest. The region with the first `split` blockers moved
    // left and the others right, for the first split that can move.
    let region: Vec<(usize, Instruction)> =
        (start..end + 1).map(|k| (k, insts[k].clone())).collect();
    let moved = |first: usize, last: usize| -> Option<Vec<(usize, Instruction)>> {
        (0..=blockers.len()).rev().find_map(|split| {
            let mut work = region.clone();
            for &b in &blockers[..split] {
                if !move_blocker(&mut work, b, first, pair, true) {
                    return None;
                }
            }
            for &b in blockers[split..].iter().rev() {
                if !move_blocker(&mut work, b, last, pair, false) {
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
        let crossed = range
            .filter(|&k| is_member(k))
            .map(|k| insts[k].clone())
            .collect();
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
    let mut best: Option<(usize, Vec<(usize, Instruction)>)> = None;
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
    let g = start
        + work
            .iter()
            .position(|(k, _)| *k == i)
            .expect("the gate itself");
    insts.splice(start..end + 1, work.into_iter().map(|(_, x)| x));
    Some(g + 1)
}

/// Move the blocker with index `blocker` in `work` (instructions tagged with
/// their index in the circuit) to just before the instruction tagged
/// `anchor`, or just after it if not `leftward`. It must commute with the
/// product of the gates on `pair` it crosses. The gates on its other qubit
/// that it depends on (or that depend on it, moving right), and in turn on
/// theirs, move with it: they share no qubit with the gates they cross, so
/// only the order of the list changes for them, not the circuit.
fn move_blocker(
    work: &mut Vec<(usize, Instruction)>,
    blocker: usize,
    anchor: usize,
    pair: [usize; 2],
    leftward: bool,
) -> bool {
    let position = |work: &[(usize, Instruction)], tag: usize| {
        work.iter()
            .position(|(k, _)| *k == tag)
            .expect("a tagged instruction")
    };
    let (p, a) = (position(work, blocker), position(work, anchor));
    let gate = work[p].1.clone();
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
        let x = &work[k].1;
        if touches(x, pair) {
            if !on_pair(x, pair) {
                return false;
            }
            crossed.push(x.clone());
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
    if !crosses(&gate, crossed, pair) {
        return false;
    }
    carried.sort_unstable();
    let group: Vec<(usize, Instruction)> = carried.iter().map(|&k| work[k].clone()).collect();
    let mut rest: Vec<(usize, Instruction)> = work
        .iter()
        .enumerate()
        .filter(|(k, _)| !carried.contains(k))
        .map(|(_, x)| x.clone())
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
fn crosses(blocker: &Instruction, crossed: Vec<Instruction>, pair: [usize; 2]) -> bool {
    let Some(shared) = blocker.qubits.iter().find(|q| pair.contains(q)).copied() else {
        return false;
    };
    let Some(other) = blocker.qubits.iter().find(|&&q| q != shared).copied() else {
        return false;
    };
    // Only the gates acting on the shared qubit can fail to commute.
    if !crossed.iter().any(|x| x.qubits.contains(&shared)) {
        return true;
    }
    let Some(s) = (Block {
        pair,
        instructions: crossed,
    })
    .matrix() else {
        return false;
    };
    let unshared = if pair[0] == shared { pair[1] } else { pair[0] };
    let qubits = [shared, unshared, other];
    let s = embed(&s, [position(&qubits, pair[0]), position(&qubits, pair[1])]);
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testing::{Rng, circuit_unitary, dense_equal_up_to_phase};

    fn gate(name: &str, qubits: &[usize], params: &[f64]) -> Instruction {
        Instruction::gate(name, qubits.to_vec(), params.to_vec())
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
}
