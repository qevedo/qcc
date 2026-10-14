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
//!    each.
//! 2. [`Timeslices`] turns a circuit into its two-qubit gates grouped into
//!    steps that can run in parallel.
//! 3. [`initial_placement`] places the qubits before the circuit starts.
//! 4. A mapper gives an assignment of qubits to cores for every timeslice, as
//!    a [`Mapping`], whose [`Mapping::communications`] is the cost:
//!    * [`hqa`] — Hungarian Qubit Assignment (Escofet et al.), the main one;
//!    * [`naive`] and [`fgp_roee`] — the baselines it is compared with.
//! 5. [`random_circuit_bounds`] gives the theoretical bounds for random
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
//! | `bounds.rs` | theoretical bounds for random circuits |

mod architecture;
mod assignment;
mod baselines;
mod bounds;
mod hqa;
pub mod hungarian;
mod lookahead;
pub mod oee;
mod timeslices;

pub use architecture::Architecture;
pub use assignment::{Assignment, Mapping, in_order, moves};
pub use baselines::{fgp_roee, naive};
pub use bounds::random_circuit_bounds;
pub use hqa::hqa;
pub use lookahead::{HORIZON, attraction};
pub use timeslices::Timeslices;

/// A mapping algorithm.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Method {
    /// Hungarian Qubit Assignment with look-ahead attraction (the default).
    Hqa,
    /// HQA without the attraction term.
    HqaWithoutLookahead,
    /// Fine-Grained Partitioning with relaxed OEE.
    FgpRoee,
    /// One qubit of each split gate moves to the other's core.
    Naive,
}

impl Method {
    pub fn parse(name: &str) -> Result<Method, String> {
        Ok(match name {
            "hqa" => Method::Hqa,
            "hqa-no-lookahead" => Method::HqaWithoutLookahead,
            "fgp-roee" => Method::FgpRoee,
            "naive" => Method::Naive,
            other => {
                return Err(format!(
                    "unknown mapping method '{other}' (hqa, hqa-no-lookahead, fgp-roee or naive)"
                ));
            }
        })
    }
}

/// Place the qubits before the circuit starts: an OEE partition of the total
/// interaction graph (how often each pair of qubits interacts), as in the
/// papers. Empty physical slots take part as qubits without interactions.
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
    part
}

/// Map `slices` onto `arch` with `method`, from [`initial_placement`].
///
/// Every slice must fit the cores (see [`Timeslices::fitted_to`]).
pub fn map(slices: &Timeslices, arch: &Architecture, method: Method) -> Result<Mapping, String> {
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
        Method::Hqa => hqa(slices, arch, &initial, true),
        Method::HqaWithoutLookahead => hqa(slices, arch, &initial, false),
        Method::FgpRoee => fgp_roee(slices, arch, &initial),
        Method::Naive => naive(slices, arch, &initial, 0x5EED),
    }?;
    debug_assert!(mapping.validate(arch, slices).is_ok());
    Ok(mapping)
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
            let arch = Architecture::new(cores, per_core).unwrap();
            let slices = Timeslices::from_pairs(qubits, &pairs)
                .unwrap()
                .fitted_to(&arch);
            for method in [
                Method::Hqa,
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
