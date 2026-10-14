//! The mappers HQA is compared with in Escofet et al.
//!
//! * [`naive`]: for each unfeasible gate, move one qubit to the other's core,
//!   swapping a random idle qubit out when that core is full (Algorithm 1 of
//!   the 2025 paper). With small cores of odd size, when that is impossible,
//!   the other qubit moves instead, or both go to a third core.
//! * [`fgp_roee`]: Fine-Grained Partitioning with relaxed OEE (Baker et al.,
//!   "Time-sliced quantum circuit partitioning for modular architectures",
//!   2020): before each slice, the qubits are re-partitioned with OEE on the
//!   look-ahead graph, where pairs that interact in the slice must end up
//!   together, starting from the previous assignment and stopping as soon as
//!   the partition is valid. With small cores of odd size rOEE can get
//!   stuck; the gates still split are then fixed as the naive mapper does.

use super::architecture::Architecture;
use super::assignment::{Assignment, Mapping, occupancy};
use super::lookahead::weights;
use super::oee;
use super::timeslices::Timeslices;
use crate::rng::Rng;

/// The naive mapper, with its random choices seeded by `seed`.
pub fn naive(
    slices: &Timeslices,
    arch: &Architecture,
    initial: &Assignment,
    seed: u64,
) -> Result<Mapping, String> {
    arch.check_fits(slices.qubits())?;
    let mut rng = Rng::new(seed);
    let mut current = initial.clone();
    let mut assignments = Vec::with_capacity(slices.len());
    for t in 0..slices.len() {
        fix_split_gates(slices, arch, t, &mut current, &mut rng)?;
        assignments.push(current.clone());
    }
    Ok(Mapping {
        architecture: arch.clone(),
        initial: initial.clone(),
        slices: assignments,
    })
}

/// Bring the qubits of every gate of slice `t` that is split across cores
/// together, the naive way: b goes to a's core, swapping a random qubit out
/// when it is full. With small cores of odd size that core can be full of
/// gates that must stay, so then the other direction is tried, then both
/// qubits into a third core.
fn fix_split_gates(
    slices: &Timeslices,
    arch: &Architecture,
    t: usize,
    current: &mut [usize],
    rng: &mut Rng,
) -> Result<(), String> {
    for &(a, b) in slices.slice(t) {
        if current[a] == current[b] {
            continue;
        }
        // A failed `bring` changes nothing, so both cores stay as they were.
        let (core_a, core_b) = (current[a], current[b]);
        let moved = bring(slices, arch, t, current, b, core_a, &[a], rng)
            || bring(slices, arch, t, current, a, core_b, &[b], rng)
            || (0..arch.cores).any(|core| {
                room(slices, arch, t, current, core, &[]) >= 2
                    && bring(slices, arch, t, current, a, core, &[], rng)
                    && bring(slices, arch, t, current, b, core, &[a], rng)
            });
        if !moved {
            return Err(format!(
                "timeslice {t}: no core can take qubits {a} and {b} together"
            ));
        }
    }
    Ok(())
}

/// Qubits of `core` that may be swapped out, other than `keep`: idle ones,
/// or else ones whose own gate in slice `t` is split across cores anyway.
fn swappable(
    slices: &Timeslices,
    t: usize,
    current: &[usize],
    core: usize,
    keep: &[usize],
) -> Vec<usize> {
    let in_core = |q: usize| current[q] == core && !keep.contains(&q);
    let idle: Vec<usize> = (0..current.len())
        .filter(|&q| in_core(q) && slices.partner(t, q).is_none())
        .collect();
    if !idle.is_empty() {
        return idle;
    }
    (0..current.len())
        .filter(|&q| in_core(q) && slices.partner(t, q).is_some_and(|p| current[p] != core))
        .collect()
}

/// How many qubits `core` could take: its free slots plus what can be swapped out.
fn room(
    slices: &Timeslices,
    arch: &Architecture,
    t: usize,
    current: &[usize],
    core: usize,
    keep: &[usize],
) -> usize {
    let free = arch.qubits_per_core - occupancy(arch, current)[core];
    free + swappable(slices, t, current, core, keep).len()
}

/// Move `q` into `core`, swapping a random swappable qubit (not in `keep`) out
/// to `q`'s old core when `core` is full. False when that is impossible.
#[allow(clippy::too_many_arguments)]
fn bring(
    slices: &Timeslices,
    arch: &Architecture,
    t: usize,
    current: &mut [usize],
    q: usize,
    core: usize,
    keep: &[usize],
    rng: &mut Rng,
) -> bool {
    if current[q] == core {
        return true;
    }
    if occupancy(arch, current)[core] == arch.qubits_per_core {
        let candidates = swappable(slices, t, current, core, keep);
        if candidates.is_empty() {
            return false;
        }
        current[candidates[rng.below(candidates.len())]] = current[q];
    }
    current[q] = core;
    true
}

/// `assignment` followed by one node per empty slot, in the cores that have them.
fn with_empty_slots(arch: &Architecture, assignment: &[usize], total: usize) -> Vec<usize> {
    let mut part = assignment.to_vec();
    let mut slots = occupancy(arch, assignment);
    for _ in assignment.len()..total {
        let core = (0..arch.cores)
            .find(|&c| slots[c] < arch.qubits_per_core)
            .expect("room left");
        slots[core] += 1;
        part.push(core);
    }
    part
}

/// Edges between qubits that interact in the current slice: heavy enough to
/// outweigh any sum of look-ahead weights (each at most 1 per qubit).
const MUST_STAY_TOGETHER: f64 = 1e4;

/// FGP-rOEE.
pub fn fgp_roee(
    slices: &Timeslices,
    arch: &Architecture,
    initial: &Assignment,
) -> Result<Mapping, String> {
    let n = slices.qubits();
    arch.check_fits(n)?;
    // Empty slots are nodes without edges, so every core is exactly full.
    let total = arch.capacity();
    let mut part = with_empty_slots(arch, initial, total);
    let mut rng = Rng::new(0x5EED);
    let mut assignments = Vec::with_capacity(slices.len());
    for t in 0..slices.len() {
        let mut w = vec![vec![0.0; total]; total];
        for (row, lookahead) in w.iter_mut().zip(weights(slices, t)) {
            row[..n].copy_from_slice(&lookahead);
        }
        for &(a, b) in slices.slice(t) {
            w[a][b] = MUST_STAY_TOGETHER;
            w[b][a] = MUST_STAY_TOGETHER;
        }
        let pairs = slices.slice(t).to_vec();
        let valid = move |p: &[usize]| pairs.iter().all(|&(a, b)| p[a] == p[b]);
        oee::partition(&w, &mut part, arch.cores, Some(&valid));
        if !valid(&part) {
            // rOEE can get stuck with small cores of odd size: fix the gates
            // still split the naive way, then refill the empty slots.
            let mut current = part[..n].to_vec();
            fix_split_gates(slices, arch, t, &mut current, &mut rng)?;
            part = with_empty_slots(arch, &current, total);
        }
        assignments.push(part[..n].to_vec());
    }
    Ok(Mapping {
        architecture: arch.clone(),
        initial: initial.clone(),
        slices: assignments,
    })
}

#[cfg(test)]
mod tests {
    use super::super::{
        hqa::{HqaOptions, hqa},
        initial_placement,
    };
    use super::*;

    #[test]
    fn baselines_are_valid_and_hqa_beats_them_on_random_circuits() {
        let mut rng = Rng::new(7);
        let arch = Architecture::new(5, 6).unwrap();
        let qubits = 30;
        let (mut naive_total, mut roee_total, mut hqa_total) = (0, 0, 0);
        for _ in 0..4 {
            let pairs: Vec<(usize, usize)> = (0..8 * qubits)
                .map(|_| {
                    let a = rng.below(qubits);
                    (a, (a + 1 + rng.below(qubits - 1)) % qubits)
                })
                .collect();
            let slices = Timeslices::from_pairs(qubits, &pairs).unwrap();
            let initial = initial_placement(&slices, &arch);
            let n = naive(&slices, &arch, &initial, 1).unwrap();
            let r = fgp_roee(&slices, &arch, &initial).unwrap();
            let h = hqa(&slices, &arch, &initial, HqaOptions::default()).unwrap();
            for m in [&n, &r, &h] {
                m.validate(&arch, &slices).unwrap();
            }
            naive_total += n.communications();
            roee_total += r.communications();
            hqa_total += h.communications();
        }
        assert!(
            hqa_total < naive_total,
            "HQA {hqa_total} vs naive {naive_total}"
        );
        assert!(
            hqa_total < roee_total,
            "HQA {hqa_total} vs FGP-rOEE {roee_total}"
        );
    }
}
