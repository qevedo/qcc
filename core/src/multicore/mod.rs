//! Mapping circuits onto multi-core (modular) quantum computers.
//!
//! A multi-core quantum computer is several small processors ("cores") linked
//! by quantum channels. A two-qubit gate needs both qubits in the same core,
//! and moving a qubit to another core (by teleportation) is slow and noisy,
//! so a mapper decides, before every step of the circuit, which core each
//! qubit is in, using as few moves (*non-local communications*) as possible.
//!
//! # How it fits together
//!
//! 1. [`Architecture`] describes the machine: how many cores, how many qubits
//!    each, and how the cores are linked (all to all, line, ring, grid or any
//!    connected set of links).
//! 2. [`Timeslices`] turns a circuit into its two-qubit gates grouped into
//!    steps that can run in parallel.
//! 3. [`initial_placement`] places the qubits before the circuit starts.
//! 4. A mapper gives an assignment of qubits to cores for every timeslice, as
//!    a [`Mapping`], whose [`Mapping::communications`] is the cost:
//!    * [`hqa`] — Hungarian Qubit Assignment (Escofet et al.), the main one;
//!    * [`naive`] and [`fgp_roee`] — the baselines it is compared with.
//! 5. [`physical_circuit`] turns a mapping into instructions on physical
//!    qubits when the cores are not all-to-all inside: SWAPs routed with
//!    SABRE inside each core and teleportations between cores
//!    ([`map_and_route`] does steps 2 to 5 at once).
//! 6. [`schedule`] works out how long that circuit takes when each link
//!    makes EPR pairs at a limited rate, on a limited number of channels.
//! 7. [`random_circuit_bounds`] gives the theoretical bounds for random
//!    circuits, to put the numbers in context.
//!
//! [`map`] runs the whole pipeline for a [`Method`].
//!
//! # Files
//!
//! | file | what it holds |
//! |---|---|
//! | `architecture.rs` | the hardware model |
//! | `timeslices.rs` | circuits as timeslices of two-qubit gates |
//! | `assignment.rs` | assignments, mappings, validity, communication count |
//! | `lookahead.rs` | how future interactions attract qubits to cores |
//! | `hungarian.rs` | the Hungarian algorithm (minimum-cost assignment) |
//! | `oee.rs` | Overall Extreme Exchange graph partitioning |
//! | `hqa.rs` | Hungarian Qubit Assignment |
//! | `baselines.rs` | the naive and FGP-rOEE mappers |
//! | `physical.rs` | routing inside cores and teleportations: a circuit that runs |
//! | `schedule.rs` | timing: how long that circuit takes on links of limited capacity |
//! | `bounds.rs` | theoretical bounds for random circuits |

mod architecture;
mod assignment;
mod baselines;
mod bounds;
mod hqa;
pub mod hungarian;
mod lookahead;
pub mod oee;
mod physical;
mod schedule;
mod timeslices;

pub use architecture::Architecture;
pub use assignment::{Assignment, Mapping, in_order, moves};
pub use baselines::{fgp_roee, naive};
pub use bounds::random_circuit_bounds;
pub use hqa::{HqaOptions, hqa};
pub use lookahead::{HORIZON, attraction};
pub use physical::{PhysicalCircuit, map_and_route, physical_circuit};
pub use schedule::{Schedule, Timing, mapping_duration, schedule};
pub use timeslices::Timeslices;

/// The weights of the delay term (in teleportations per teleportation of
/// delay) that [`Method::HqaLatency`] tries; 0 is plain HQA.
pub const LATENCY_WEIGHTS: [f64; 4] = [0.0, 0.25, 0.5, 1.0];

/// A mapping algorithm.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Method {
    /// Hungarian Qubit Assignment with look-ahead attraction (the default).
    Hqa,
    /// HQA that also weighs how much later each placement lets its gate start
    /// (not in the papers), with the weight from [`LATENCY_WEIGHTS`] whose
    /// mapping takes least time under the [`Timing`] given to [`map_timed`].
    HqaLatency,
    /// HQA without the attraction term.
    HqaWithoutLookahead,
    /// HQA that counts every move as one teleportation, ignoring how many
    /// links it crosses (the papers' cost, for comparison on sparse links).
    HqaUniformCost,
    /// Fine-Grained Partitioning with relaxed OEE.
    FgpRoee,
    /// One qubit of each split gate moves to the other's core.
    Naive,
}

impl Method {
    pub fn parse(name: &str) -> Result<Method, String> {
        Ok(match name {
            "hqa" => Method::Hqa,
            "hqa-latency" => Method::HqaLatency,
            "hqa-no-lookahead" => Method::HqaWithoutLookahead,
            "hqa-uniform-cost" => Method::HqaUniformCost,
            "fgp-roee" => Method::FgpRoee,
            "naive" => Method::Naive,
            other => {
                return Err(format!(
                    "unknown mapping method '{other}' (hqa, hqa-latency, hqa-no-lookahead, hqa-uniform-cost, fgp-roee or naive)"
                ));
            }
        })
    }
}

/// Place the qubits before the circuit starts: an OEE partition of the total
/// interaction graph (how often each pair of qubits interacts), as in the
/// papers. Empty physical slots take part as qubits without interactions.
/// When the cores are not all linked, the groups that interact most go to
/// nearby cores.
pub fn initial_placement(slices: &Timeslices, arch: &Architecture) -> Assignment {
    let n = slices.qubits();
    let total = arch.capacity().max(n);
    let counts = slices.interaction_counts();
    let mut w = vec![vec![0.0; total]; total];
    for (row, counts) in w.iter_mut().zip(counts) {
        row[..n].copy_from_slice(&counts);
    }
    let mut part = in_order(arch, total);
    oee::partition(&w, &mut part, arch.cores, None);
    part.truncate(n);
    if arch.is_all_to_all() {
        return part;
    }
    // On sparse links, put the groups that interact most on nearby cores.
    let order = place_groups(&slices.interaction_counts(), &part, arch);
    part.into_iter().map(|group| order[group]).collect()
}

/// For OEE groups `part` on a sparse `arch`, the core of each group that
/// lowers Σ (interactions between two groups × links between their cores):
/// start in order and swap pairs of groups while that helps.
fn place_groups(counts: &[Vec<f64>], part: &[usize], arch: &Architecture) -> Vec<usize> {
    let k = arch.cores;
    let mut traffic = vec![vec![0.0; k]; k];
    for (a, row) in counts.iter().enumerate() {
        for (b, &w) in row.iter().enumerate() {
            traffic[part[a]][part[b]] += w;
        }
    }
    let mut core: Vec<usize> = (0..k).collect();
    let cost = |core: &[usize]| -> f64 {
        (0..k)
            .flat_map(|g| (0..k).map(move |h| (g, h)))
            .map(|(g, h)| traffic[g][h] * arch.hops(core[g], core[h]) as f64)
            .sum()
    };
    let mut best = cost(&core);
    loop {
        let mut improved = false;
        for g in 0..k {
            for h in g + 1..k {
                core.swap(g, h);
                let c = cost(&core);
                if c + 1e-9 < best {
                    best = c;
                    improved = true;
                } else {
                    core.swap(g, h);
                }
            }
        }
        if !improved {
            return core;
        }
    }
}

/// Map `slices` onto `arch` with `method`, from [`initial_placement`].
///
/// Every slice must fit the cores (see [`Timeslices::fitted_to`]).
/// [`Method::HqaLatency`] assumes the default [`Timing`]; see [`map_timed`].
pub fn map(slices: &Timeslices, arch: &Architecture, method: Method) -> Result<Mapping, String> {
    map_timed(slices, arch, method, &Timing::default())
}

/// [`map`], with the durations and link capacity that
/// [`Method::HqaLatency`] plans for (the other methods ignore them).
pub fn map_timed(
    slices: &Timeslices,
    arch: &Architecture,
    method: Method,
    timing: &Timing,
) -> Result<Mapping, String> {
    arch.check_fits(slices.qubits())?;
    if !slices.fits(arch) {
        return Err(format!(
            "a timeslice has more gates than {} cores of {} qubits can run at once: \
             split it with Timeslices::fitted_to",
            arch.cores, arch.qubits_per_core
        ));
    }
    let initial = initial_placement(slices, arch);
    let mapping = match method {
        Method::Hqa => hqa(slices, arch, &initial, HqaOptions::default()),
        Method::HqaLatency => hqa_latency(slices, arch, &initial, timing),
        Method::HqaWithoutLookahead => hqa(
            slices,
            arch,
            &initial,
            HqaOptions {
                lookahead: false,
                ..HqaOptions::default()
            },
        ),
        Method::HqaUniformCost => hqa(
            slices,
            arch,
            &initial,
            HqaOptions {
                hops: false,
                ..HqaOptions::default()
            },
        ),
        Method::FgpRoee => fgp_roee(slices, arch, &initial),
        Method::Naive => naive(slices, arch, &initial, 0x5EED),
    }?;
    debug_assert!(mapping.validate(arch, slices).is_ok());
    Ok(mapping)
}

/// HQA with each of [`LATENCY_WEIGHTS`]: the mapping that takes least time
/// under `timing` (then the one with fewest communications). Weighing delay
/// shortens circuits whose links have pairs to spare, but its extra
/// teleportations lengthen those that keep the links busy; trying the
/// weights, plain HQA included, finds which case a circuit is in.
fn hqa_latency(
    slices: &Timeslices,
    arch: &Architecture,
    initial: &Assignment,
    timing: &Timing,
) -> Result<Mapping, String> {
    let candidates = latency_candidates(slices, arch, initial, timing)?;
    let key = |m: &Mapping| {
        (
            schedule::mapping_duration(m, slices, timing),
            m.communications(),
        )
    };
    Ok(fastest(candidates, key))
}

/// HQA's mappings with each of [`LATENCY_WEIGHTS`].
fn latency_candidates(
    slices: &Timeslices,
    arch: &Architecture,
    initial: &Assignment,
    timing: &Timing,
) -> Result<Vec<Mapping>, String> {
    LATENCY_WEIGHTS
        .iter()
        .map(|&latency| {
            let options = HqaOptions {
                latency,
                teleport_time: timing.teleport / timing.two_qubit,
                ..HqaOptions::default()
            };
            hqa(slices, arch, initial, options)
        })
        .collect()
}

/// The item with the smallest `(duration, count)`, the first on ties.
fn fastest<T>(items: Vec<T>, key: impl Fn(&T) -> (f64, usize)) -> T {
    let mut best: Option<((f64, usize), T)> = None;
    for item in items {
        let k = key(&item);
        if best
            .as_ref()
            .is_none_or(|(b, _)| k.0 < b.0 - 1e-9 || (k.0 <= b.0 + 1e-9 && k.1 < b.1))
        {
            best = Some((k, item));
        }
    }
    best.expect("at least one item").1
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rng::Rng;

    /// Random architectures (even and odd core sizes, full and with room) and
    /// random circuits: every method must give a valid mapping.
    #[test]
    fn every_method_maps_any_circuit_validly() {
        let mut rng = Rng::new(99);
        for _ in 0..300 {
            let cores = 2 + rng.below(6);
            let per_core = 2 + rng.below(7);
            let qubits = 2 + rng.below(cores * per_core - 1);
            let pairs: Vec<(usize, usize)> = (0..1 + rng.below(8 * qubits))
                .map(|_| {
                    let a = rng.below(qubits);
                    (a, (a + 1 + rng.below(qubits - 1)) % qubits)
                })
                .collect();
            let arch = match rng.below(4) {
                0 => Architecture::new(cores, per_core),
                1 => Architecture::line(cores, per_core),
                2 => Architecture::ring(cores, per_core),
                _ => Architecture::grid(2, cores.div_ceil(2), per_core),
            }
            .unwrap();
            let slices = Timeslices::from_pairs(qubits, &pairs)
                .unwrap()
                .fitted_to(&arch);
            for method in [
                Method::Hqa,
                Method::HqaLatency,
                Method::HqaWithoutLookahead,
                Method::FgpRoee,
                Method::Naive,
            ] {
                let mapping = map(&slices, &arch, method)
                    .unwrap_or_else(|e| panic!("{method:?} {cores}×{per_core}: {e}"));
                mapping.validate(&arch, &slices).unwrap();
            }
        }
    }

    #[test]
    fn oversized_slices_are_refused_until_fitted() {
        let arch = Architecture::new(2, 3).unwrap();
        let slices = Timeslices::from_pairs(6, &[(0, 3), (1, 4), (2, 5)]).unwrap();
        assert!(
            map(&slices, &arch, Method::Hqa)
                .unwrap_err()
                .contains("fitted_to")
        );
        assert!(map(&slices.fitted_to(&arch), &arch, Method::Hqa).is_ok());
    }
}

#[cfg(test)]
mod topology_tests {
    use super::*;
    use crate::rng::Rng;

    fn random_layers(rng: &mut Rng, qubits: usize, layers: usize) -> Vec<(usize, usize)> {
        let mut pairs = Vec::new();
        for _ in 0..layers {
            let mut order: Vec<usize> = (0..qubits).collect();
            for i in (1..qubits).rev() {
                order.swap(i, rng.below(i + 1));
            }
            pairs.extend(order.chunks(2).map(|c| (c[0], c[1])));
        }
        pairs
    }

    /// On a line of cores, counting links in the cost (HQA) beats counting
    /// every move as one (the papers' cost).
    #[test]
    fn counting_links_helps_on_a_line() {
        let mut rng = Rng::new(4);
        let arch = Architecture::line(8, 6).unwrap();
        let (mut hops, mut uniform) = (0, 0);
        for _ in 0..4 {
            let slices = Timeslices::from_pairs(48, &random_layers(&mut rng, 48, 40)).unwrap();
            hops += map(&slices, &arch, Method::Hqa).unwrap().communications();
            uniform += map(&slices, &arch, Method::HqaUniformCost)
                .unwrap()
                .communications();
        }
        assert!(hops < uniform, "with links {hops}, uniform {uniform}");
    }

    #[test]
    fn all_to_all_communications_are_moves() {
        let mut rng = Rng::new(5);
        let arch = Architecture::new(4, 6).unwrap();
        let slices = Timeslices::from_pairs(24, &random_layers(&mut rng, 24, 10)).unwrap();
        let mapping = map(&slices, &arch, Method::Hqa).unwrap();
        assert_eq!(mapping.communications(), mapping.moved_qubits());
        // And HQA's cost there is the papers' 1 or 2.
        let uniform = map(&slices, &arch, Method::HqaUniformCost).unwrap();
        assert_eq!(uniform.communications(), mapping.communications());
    }
}
