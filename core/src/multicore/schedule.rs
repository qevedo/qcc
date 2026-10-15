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
use super::assignment::Mapping;
use super::timeslices::Timeslices;
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
    let mut links = Links::new(arch, timing);
    let mut free = vec![0.0f64; arch.capacity()];
    let mut start = Vec::with_capacity(instructions.len());
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
            let t = links.teleport(arch, timing, from / k, to / k, begin);
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
        epr_wait: links.epr_wait,
        link_hops: arch.links().iter().copied().zip(links.hops).collect(),
    })
}

/// How long a mapping would take under `timing`, before routing inside the
/// cores: each two-qubit gate takes `two_qubit`, each move a teleportation
/// per link, sharing the links' channels as in [`schedule`], and everything
/// starts as soon as its qubits are free. Mappers use it to compare
/// mappings by time.
pub fn mapping_duration(mapping: &Mapping, slices: &Timeslices, timing: &Timing) -> f64 {
    let arch = &mapping.architecture;
    let mut links = Links::new(arch, timing);
    let mut free = vec![0.0f64; mapping.initial.len()];
    let mut previous = &mapping.initial;
    let mut duration = 0.0f64;
    for (t, assignment) in mapping.slices.iter().enumerate() {
        for q in 0..free.len() {
            if previous[q] != assignment[q] {
                free[q] = links.teleport(arch, timing, previous[q], assignment[q], free[q]);
            }
        }
        for &(a, b) in slices.slice(t) {
            let end = free[a].max(free[b]) + timing.two_qubit;
            free[a] = end;
            free[b] = end;
            duration = duration.max(end);
        }
        previous = assignment;
    }
    free.iter().copied().fold(duration, f64::max)
}

/// The EPR pairs of every link.
///
/// A channel holds one pair and makes the next in `epr_pair` once it is used,
/// so two uses of a channel are at least `epr_pair` apart. Uses are booked in
/// the order operations are scheduled, which is not always the order in time,
/// so a use can go in a gap between earlier bookings.
struct Links {
    /// `uses[l][c]`: when channel `c` of link `l` hands out a pair, sorted.
    uses: Vec<Vec<Vec<f64>>>,
    hops: Vec<usize>,
    epr_wait: f64,
}

impl Links {
    fn new(arch: &Architecture, timing: &Timing) -> Links {
        Links {
            uses: vec![vec![Vec::new(); timing.link_capacity.max(1)]; arch.links().len()],
            hops: vec![0; arch.links().len()],
            epr_wait: 0.0,
        }
    }

    /// Teleport a qubit, free at `t`, from core `a` to core `b`, one hop per
    /// link, each with the first pair ready on it. Returns when it arrives.
    fn teleport(
        &mut self,
        arch: &Architecture,
        timing: &Timing,
        a: usize,
        b: usize,
        t: f64,
    ) -> f64 {
        let mut t = t;
        for link in path(arch, a, b) {
            let l = arch
                .links()
                .binary_search(&link)
                .expect("a path follows links");
            let channels = &mut self.uses[l];
            let (c, hop_start) = channels
                .iter()
                .map(|uses| earliest_use(uses, t, timing.epr_pair))
                .enumerate()
                .min_by(|a, b| a.1.total_cmp(&b.1))
                .expect("at least one channel");
            let uses = &mut channels[c];
            let at = uses.partition_point(|&x| x <= hop_start);
            uses.insert(at, hop_start);
            self.epr_wait += hop_start - t;
            self.hops[l] += 1;
            t = hop_start + timing.teleport;
        }
        t
    }
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

/// The earliest time from `t` when a channel with these sorted `uses` can hand
/// out a pair: at least `gap` away from every other use.
fn earliest_use(uses: &[f64], t: f64, gap: f64) -> f64 {
    let mut h = t;
    for &x in &uses[uses.partition_point(|&x| x <= t - gap)..] {
        if h + gap <= x {
            break;
        }
        h = h.max(x + gap);
    }
    h
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
    fn mapping_durations_count_gates_and_moves() {
        let arch = Architecture::new(2, 3).unwrap();
        let slices = Timeslices::from_pairs(4, &[(0, 1), (1, 2), (2, 3)]).unwrap();
        // q1 moves to q2's core for slice 1 and stays.
        let mapping = Mapping {
            architecture: arch.clone(),
            initial: vec![0, 0, 1, 1],
            slices: vec![vec![0, 0, 1, 1], vec![0, 1, 1, 1], vec![0, 1, 1, 1]],
        };
        mapping.validate(&arch, &slices).unwrap();
        // Gate (1), q1's teleportation (5), then two gates in a row (2).
        let d = mapping_duration(&mapping, &slices, &Timing::default());
        assert!((d - 8.0).abs() < 1e-12, "{d}");
    }

    #[test]
    fn teleportations_use_gaps_between_earlier_bookings() {
        let arch = Architecture::new(2, 4).unwrap();
        // The first teleportation waits for its qubit until 10; the second,
        // listed later but ready at 5, does not wait behind it.
        let circuit = vec![
            Instruction::gate("h", vec![0], vec![]),
            teleport(0, 4),
            Instruction::gate("measure", vec![1], vec![]),
            teleport(1, 5),
        ];
        let timing = Timing {
            one_qubit: 10.0,
            epr_pair: 0.0,
            ..Timing::default()
        };
        let s = schedule(&circuit, &arch, &timing).unwrap();
        assert_eq!(s.start[3], 5.0);
        assert_eq!(s.epr_wait, 0.0);
        assert_eq!(earliest_use(&[10.0], 0.0, 10.0), 0.0);
        assert_eq!(earliest_use(&[10.0], 1.0, 10.0), 20.0);
        assert_eq!(earliest_use(&[0.0, 30.0], 5.0, 10.0), 10.0);
        assert_eq!(earliest_use(&[0.0, 15.0], 5.0, 10.0), 25.0);
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
