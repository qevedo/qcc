//! How long a multi-core circuit takes, when the links between cores can only
//! teleport so many qubits at a time.
//!
//! Counting teleportations treats every one alike. On hardware, a
//! teleportation consumes an EPR pair shared across a link, and a link
//! produces pairs at a limited rate, on a limited number of channels (pairs
//! of communication qubits) at once. When many qubits cross the same link
//! together, they wait for pairs. This module schedules a
//! [`PhysicalCircuit`](super::PhysicalCircuit) under such a model and reports
//! how long it takes and how long teleportations waited.
//!
//! The model ([`Timing`]):
//!
//! - every operation takes a fixed time by kind: single-qubit gates,
//!   two-qubit gates, SWAPs, measurements;
//! - each link has `link_capacity` channels; a channel holds one EPR pair at a
//!   time and makes a new one, in `epr_pair`, as soon as its pair is used, so
//!   a channel that has been idle has a pair ready;
//! - one hop of a teleportation uses one pair of a link on the way and then
//!   takes `teleport` (Bell measurement, classical message, correction);
//!   moving between cores that are not linked takes one hop per link, one
//!   after the other, along a shortest path.
//!
//! Operations start as soon as their qubits are free (as soon as possible, in
//! the order of the circuit), and a hop as soon as a pair is ready on one of
//! the link's channels. As in [`PhysicalCircuit`](super::PhysicalCircuit),
//! a run of consecutive teleportations happens at once: a qubit arriving in a
//! slot does not wait for the qubit leaving it.

use super::architecture::Architecture;
use crate::ir::Instruction;

/// Durations and link resources, in any one time unit.
#[derive(Clone, Copy, Debug)]
pub struct Timing {
    pub one_qubit: f64,
    pub two_qubit: f64,
    pub swap: f64,
    pub measure: f64,
    /// Time for a channel to make one EPR pair.
    pub epr_pair: f64,
    /// Time of one hop of a teleportation once its pair is ready.
    pub teleport: f64,
    /// Channels per link: EPR pairs it can hold, and so teleportations it can
    /// start, at once.
    pub link_capacity: usize,
}

impl Default for Timing {
    /// Illustrative values in units of a two-qubit gate: a SWAP is three, a
    /// measurement five, an EPR pair ten and the rest of a teleportation five,
    /// with one channel per link. Set them for the hardware at hand.
    fn default() -> Timing {
        Timing {
            one_qubit: 0.1,
            two_qubit: 1.0,
            swap: 3.0,
            measure: 5.0,
            epr_pair: 10.0,
            teleport: 5.0,
            link_capacity: 1,
        }
    }
}

/// When everything happens.
#[derive(Clone, Debug)]
pub struct Schedule {
    /// `start[i]`: when instruction `i` starts.
    pub start: Vec<f64>,
    /// When the last instruction ends.
    pub duration: f64,
    /// Total time teleportation hops waited for an EPR pair, beyond waiting
    /// for their qubit.
    pub epr_wait: f64,
    /// The links, as in [`Architecture::links`], and how many hops used each.
    pub link_hops: Vec<((usize, usize), usize)>,
}

/// Schedule `instructions`, on the physical qubits of `arch` (as a
/// [`PhysicalCircuit`](super::PhysicalCircuit) numbers them), under `timing`.
pub fn schedule(
    instructions: &[Instruction],
    arch: &Architecture,
    timing: &Timing,
) -> Result<Schedule, String> {
    if timing.link_capacity == 0 {
        return Err("each link needs at least one channel".into());
    }
    let k = arch.qubits_per_core;
    let links = arch.links();
    // `pair_ready[l][c]`: when channel `c` of link `l` has a pair.
    let mut pair_ready = vec![vec![0.0f64; timing.link_capacity]; links.len()];
    let mut hops = vec![0usize; links.len()];
    let mut free = vec![0.0f64; arch.capacity()];
    let mut start = Vec::with_capacity(instructions.len());
    let mut epr_wait = 0.0;
    let mut duration = 0.0f64;
    if let Some(inst) = instructions
        .iter()
        .find(|i| i.qubits.iter().any(|&q| q >= arch.capacity()))
    {
        return Err(format!(
            "'{}' uses a qubit beyond the {} of the cores",
            inst.name,
            arch.capacity()
        ));
    }
    let mut i = 0;
    while i < instructions.len() {
        if instructions[i].name != "teleport" {
            let inst = &instructions[i];
            let begin = inst.qubits.iter().map(|&q| free[q]).fold(0.0, f64::max);
            let end = begin + duration_of(inst, timing);
            for &q in &inst.qubits {
                free[q] = end;
            }
            start.push(begin);
            duration = duration.max(end);
            i += 1;
            continue;
        }
        // A run of teleportations happens at once: each waits for its qubit
        // and for the last operation on its destination, not for the qubit
        // leaving the destination in the same run (arrivals land in the
        // core's communication qubits first).
        let run = instructions[i..]
            .iter()
            .take_while(|inst| inst.name == "teleport")
            .count();
        let mut ends = Vec::with_capacity(run);
        for inst in &instructions[i..i + run] {
            let (from, to) = (inst.qubits[0], inst.qubits[1]);
            let begin = free[from].max(free[to]);
            let mut t = begin;
            for link in path(arch, from / k, to / k) {
                let l = links.binary_search(&link).expect("a path follows links");
                let channel = &mut pair_ready[l];
                let c = (0..channel.len())
                    .min_by(|&a, &b| channel[a].total_cmp(&channel[b]))
                    .expect("at least one channel");
                let hop_start = t.max(channel[c]);
                epr_wait += hop_start - t;
                channel[c] = hop_start + timing.epr_pair;
                hops[l] += 1;
                t = hop_start + timing.teleport;
            }
            start.push(begin);
            ends.push(t);
            duration = duration.max(t);
        }
        // Sources are free once their qubit has left; destinations once it
        // has arrived.
        for (inst, &begin) in instructions[i..i + run].iter().zip(&start[i..]) {
            free[inst.qubits[0]] = begin;
        }
        for (inst, &end) in instructions[i..i + run].iter().zip(&ends) {
            free[inst.qubits[1]] = end;
        }
        i += run;
    }
    Ok(Schedule {
        start,
        duration,
        epr_wait,
        link_hops: links.iter().copied().zip(hops).collect(),
    })
}

fn duration_of(inst: &Instruction, timing: &Timing) -> f64 {
    match (inst.name.as_str(), inst.qubits.len()) {
        ("barrier", _) => 0.0,
        ("measure" | "reset", _) => timing.measure,
        ("swap", 2) => timing.swap,
        (_, 0 | 1) => timing.one_qubit,
        _ => timing.two_qubit,
    }
}

/// The links on a shortest path from core `a` to core `b`, each with the
/// smaller core first.
fn path(arch: &Architecture, a: usize, b: usize) -> Vec<(usize, usize)> {
    let mut out = Vec::new();
    let mut core = a;
    while core != b {
        let next = (0..arch.cores)
            .find(|&n| arch.hops(core, n) == 1 && arch.hops(n, b) + 1 == arch.hops(core, b))
            .expect("the links connect the cores");
        out.push((core.min(next), core.max(next)));
        core = next;
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn teleport(from: usize, to: usize) -> Instruction {
        Instruction::gate("teleport", vec![from, to], vec![])
    }

    #[test]
    fn gates_run_in_parallel_and_wait_for_their_qubits() {
        let arch = Architecture::new(2, 2).unwrap();
        let timing = Timing::default();
        let circuit = vec![
            Instruction::gate("h", vec![0], vec![]),
            Instruction::gate("cx", vec![2, 3], vec![]),
            Instruction::gate("cx", vec![0, 1], vec![]),
            Instruction::gate("swap", vec![0, 1], vec![]),
        ];
        let s = schedule(&circuit, &arch, &timing).unwrap();
        assert_eq!(s.start, vec![0.0, 0.0, 0.1, 1.1]);
        assert!((s.duration - 4.1).abs() < 1e-12);
    }

    #[test]
    fn teleportations_share_the_channels_of_a_link() {
        let arch = Architecture::new(2, 4).unwrap();
        let timing = Timing::default();
        // Three qubits cross the one link at once.
        let circuit = vec![teleport(0, 4), teleport(1, 5), teleport(2, 6)];
        let one = schedule(&circuit, &arch, &timing).unwrap();
        // The first finds a pair ready; the next ones wait for new pairs.
        assert_eq!(one.start, vec![0.0, 0.0, 0.0]);
        assert!((one.duration - 25.0).abs() < 1e-12);
        assert!((one.epr_wait - 30.0).abs() < 1e-12);
        assert_eq!(one.link_hops, vec![((0, 1), 3)]);
        let three = schedule(
            &circuit,
            &arch,
            &Timing {
                link_capacity: 3,
                ..timing
            },
        )
        .unwrap();
        assert!((three.duration - 5.0).abs() < 1e-12);
        assert_eq!(three.epr_wait, 0.0);
    }

    #[test]
    fn unlinked_cores_take_one_hop_per_link() {
        let arch = Architecture::line(3, 2).unwrap();
        let s = schedule(&[teleport(0, 4)], &arch, &Timing::default()).unwrap();
        assert!((s.duration - 10.0).abs() < 1e-12);
        assert_eq!(s.link_hops, vec![((0, 1), 1), ((1, 2), 1)]);
        assert!(schedule(&[teleport(0, 9)], &arch, &Timing::default()).is_err());
    }

    #[test]
    fn qubits_swapping_cores_do_not_wait_for_each_other() {
        let arch = Architecture::new(2, 2).unwrap();
        let timing = Timing {
            link_capacity: 2,
            ..Timing::default()
        };
        // Full cores exchange a qubit, then a gate on the arrival.
        let circuit = vec![
            Instruction::gate("x", vec![2], vec![]),
            teleport(0, 2),
            teleport(2, 0),
            Instruction::gate("cx", vec![2, 3], vec![]),
        ];
        let s = schedule(&circuit, &arch, &timing).unwrap();
        assert_eq!(s.start, vec![0.0, 0.1, 0.1, 5.1]);
    }
}
