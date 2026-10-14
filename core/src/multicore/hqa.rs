//! Hungarian Qubit Assignment (HQA), from
//! P. Escofet, A. Ovide, C. G. Almudéver, E. Alarcón and S. Abadal,
//! "Hungarian Qubit Assignment for Optimized Mapping of Quantum Circuits on
//! Multi-Core Architectures" (arXiv:2309.12182), extended in "Revisiting the
//! Mapping of Quantum Circuits: Entering the Multi-Core Era" (ACM TQC 2025,
//! arXiv:2403.17205).
//!
//! Before each timeslice, the gates whose qubits are in different cores
//! ("unfeasible" gates) are collected and their qubits taken out of their
//! cores. The gates, not the qubits, are then assigned to cores: a cost matrix
//! (gates × cores) goes to the Hungarian algorithm, which gives each core at
//! most one gate; both qubits of a gate go into the free space of its core.
//! This repeats, with the free space updated, until every gate is placed.
//!
//! The cost of putting gate (qA, qB) in core c is
//!
//! ```text
//! ∞                if c has fewer than two free slots
//! 1 − attr(op, c)  if qA or qB was in c before the slice (one teleportation)
//! 2 − attr(op, c)  otherwise (two teleportations)
//! ```
//!
//! When the cores are not all linked to each other, the papers note that the
//! cost "can be adapted". Here both terms count teleportations: the links
//! qA and qB cross to reach c now, plus, for the look-ahead, the links each
//! would cross to meet its future partners from c, weighted by `2^−(m−t)`
//! (`lookahead::future_hops_to_cores`). With every core linked to every other
//! this differs from the formula above by a constant per gate, so the
//! Hungarian algorithm makes the same choices.
//!
//! where attr(op, c) is the mean attraction of qA and qB to c (see
//! [`lookahead`](super::lookahead)), so gates go where their qubits will be
//! needed next.
//!
//! Choices the papers leave open, made here:
//! * Free space must come in pairs. When the cores' free slots cannot hold the
//!   gates two by two, the papers pair up cores with an odd number of free
//!   slots with an auxiliary gate between a non-interacting qubit of each;
//!   that qubit is the one least attracted to its core. With spare capacity
//!   the auxiliary gate is only added when it is needed. When only one of the
//!   two cores has an idle qubit (small cores of odd size), that qubit moves
//!   into the other core instead, which also makes both even. When no odd
//!   core has an idle qubit, every idle qubit is taken out, the gates are
//!   placed, and the idle qubits go back to their cores where there is room.
//! * The papers' cores have an even size. With an odd size, a timeslice can
//!   hold more gates than the cores can run at once; [`Timeslices::fitted_to`]
//!   splits such slices first.
//! * Ties in the Hungarian algorithm go to the lower-numbered gate and core.
//! * Attraction looks [`HORIZON`](super::lookahead::HORIZON) slices ahead,
//!   past which its terms are below rounding.

use super::architecture::Architecture;
use super::assignment::{Assignment, Mapping, occupancy};
use super::hungarian::{self, FORBIDDEN};
use super::lookahead::{attraction, attraction_to_cores, future_hops_to_cores};
use super::timeslices::Timeslices;

/// On sparse links, moving both qubits of a gate to a core between them
/// costs as many teleportations as moving one to the other, but disturbs two
/// placements instead of one. This small cost per moved qubit breaks such
/// ties towards moving fewer qubits, as the papers' "1 if a qubit is already
/// there" does with all cores linked to all. (It is a tie-breaker: on the
/// benchmarks, values from 0.001 to 0.05 all give within 1% of each other.)
const MOVE_TIE_BREAK: f64 = 0.01;

/// Which terms the HQA cost uses.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct HqaOptions {
    /// Attraction to where qubits are needed next (the papers' default).
    pub lookahead: bool,
    /// Count the links a qubit crosses to reach a core. Without it, every
    /// move costs one, as with all cores linked to all.
    pub hops: bool,
}

impl Default for HqaOptions {
    fn default() -> HqaOptions {
        HqaOptions {
            lookahead: true,
            hops: true,
        }
    }
}

/// Map the circuit's timeslices onto the cores with HQA, starting from
/// `initial`.
pub fn hqa(
    slices: &Timeslices,
    arch: &Architecture,
    initial: &Assignment,
    options: HqaOptions,
) -> Result<Mapping, String> {
    arch.check_fits(slices.qubits())?;
    let mut current = initial.clone();
    let mut assignments = Vec::with_capacity(slices.len());
    for t in 0..slices.len() {
        let unfeasible: Vec<(usize, usize)> = slices
            .slice(t)
            .iter()
            .filter(|&&(a, b)| current[a] != current[b])
            .copied()
            .collect();
        if !unfeasible.is_empty() {
            current = place(slices, arch, t, &current, unfeasible, options)?;
        }
        assignments.push(current.clone());
    }
    Ok(Mapping {
        architecture: arch.clone(),
        initial: initial.clone(),
        slices: assignments,
    })
}

/// The assignment for slice `t`: `previous` with the qubits of the unfeasible
/// gates re-placed.
fn place(
    slices: &Timeslices,
    arch: &Architecture,
    t: usize,
    previous: &[usize],
    mut gates: Vec<(usize, usize)>,
    options: HqaOptions,
) -> Result<Assignment, String> {
    let mut core_of: Vec<Option<usize>> = previous.iter().map(|&c| Some(c)).collect();
    let mut free: Vec<usize> = occupancy(arch, previous)
        .iter()
        .map(|&used| arch.qubits_per_core - used)
        .collect();
    for &(a, b) in &gates {
        for q in [a, b] {
            free[previous[q]] += 1;
            core_of[q] = None;
        }
    }
    // Idle qubits taken out to make room, placed again after the gates.
    let mut singles = Vec::new();
    let real_gates = gates.len();
    if !balance_free_space(slices, t, &mut gates, &mut core_of, &mut free) {
        // Auxiliary gates added so far go too: their qubits become singles.
        for (a, b) in gates.drain(real_gates..) {
            singles.extend([a, b]);
        }
        // Last resort (small cores of odd size): take every idle qubit out.
        // Then each core's free space is its size minus its feasible gates,
        // which holds the unfeasible gates two by two because the slice fits
        // the cores (`Timeslices::fitted_to`).
        for q in 0..core_of.len() {
            if let (Some(c), None) = (core_of[q], slices.partner(t, q)) {
                core_of[q] = None;
                free[c] += 1;
                singles.push(q);
            }
        }
    }
    place_gates(
        slices,
        arch,
        t,
        previous,
        gates,
        &mut core_of,
        &mut free,
        options,
    )?;
    for q in singles {
        // Back to its core when there is room, else where it is drawn most.
        let mut attr = vec![0.0; arch.cores];
        attraction_to_cores(slices, t, q, &core_of, &mut attr);
        let core = if free[previous[q]] > 0 {
            previous[q]
        } else {
            (0..arch.cores)
                .filter(|&c| free[c] > 0)
                .max_by(|&x, &y| attr[x].total_cmp(&attr[y]).then(y.cmp(&x)))
                .ok_or_else(|| format!("timeslice {t}: no room left for qubit {q}"))?
        };
        core_of[q] = Some(core);
        free[core] -= 1;
    }
    Ok(core_of
        .into_iter()
        .map(|c| c.expect("every qubit is placed"))
        .collect())
}

/// Assign `gates` to cores with the Hungarian algorithm, one gate per core per
/// round, until all are placed.
#[allow(clippy::too_many_arguments)]
fn place_gates(
    slices: &Timeslices,
    arch: &Architecture,
    t: usize,
    previous: &[usize],
    mut gates: Vec<(usize, usize)>,
    core_of: &mut [Option<usize>],
    free: &mut [usize],
    options: HqaOptions,
) -> Result<(), String> {
    let mut attr_a = vec![0.0; arch.cores];
    let mut attr_b = vec![0.0; arch.cores];
    let sparse = options.hops && !arch.is_all_to_all();
    while !gates.is_empty() {
        let mut cost = Vec::with_capacity(gates.len());
        for &(a, b) in &gates {
            // The look-ahead term to subtract: the papers' attraction, or, on
            // sparse links, minus the teleportations still to come.
            if options.lookahead && sparse {
                future_hops_to_cores(slices, t, a, core_of, arch, &mut attr_a);
                future_hops_to_cores(slices, t, b, core_of, arch, &mut attr_b);
                attr_a
                    .iter_mut()
                    .chain(attr_b.iter_mut())
                    .for_each(|x| *x = -*x);
            } else if options.lookahead {
                attraction_to_cores(slices, t, a, core_of, &mut attr_a);
                attraction_to_cores(slices, t, b, core_of, &mut attr_b);
            }
            let row: Vec<f64> = (0..arch.cores)
                .map(|c| {
                    if free[c] < 2 {
                        return FORBIDDEN;
                    }
                    // Teleportations to bring both qubits to c: 1 or 2 with
                    // all cores linked to all, as in the papers.
                    let hops = |q: usize| {
                        if previous[q] == c {
                            0
                        } else if options.hops {
                            arch.hops(previous[q], c)
                        } else {
                            1
                        }
                    };
                    let moved = usize::from(previous[a] != c) + usize::from(previous[b] != c);
                    let base = (hops(a) + hops(b)) as f64
                        + if sparse {
                            MOVE_TIE_BREAK * moved as f64
                        } else {
                            0.0
                        };
                    base - (attr_a[c] + attr_b[c]) / 2.0
                })
                .collect();
            cost.push(row);
        }
        let chosen = hungarian::assign(&cost);
        let mut remaining = Vec::new();
        for (&(a, b), core) in gates.iter().zip(chosen) {
            match core {
                Some(c) => {
                    core_of[a] = Some(c);
                    core_of[b] = Some(c);
                    free[c] -= 2;
                }
                None => remaining.push((a, b)),
            }
        }
        if remaining.len() == gates.len() {
            return Err(format!(
                "timeslice {t}: no core has room for the remaining gates"
            ));
        }
        gates = remaining;
    }
    Ok(())
}

/// Make the free space hold the gates two by two, as the papers do: while it
/// cannot, take an idle qubit out of each of two cores with an odd number of
/// free slots and add an auxiliary gate between them. When only one of the
/// two cores has an idle qubit, move it into the other instead. Returns false
/// when neither works (no odd core has an idle qubit).
fn balance_free_space(
    slices: &Timeslices,
    t: usize,
    gates: &mut Vec<(usize, usize)>,
    core_of: &mut [Option<usize>],
    free: &mut [usize],
) -> bool {
    while free.iter().map(|f| f / 2).sum::<usize>() < gates.len() {
        let odd: Vec<usize> = (0..free.len()).filter(|&c| free[c] % 2 == 1).collect();
        let idle: Vec<Option<usize>> = odd
            .iter()
            .map(|&c| least_attached_idle(slices, t, c, core_of))
            .collect();
        let with_idle: Vec<usize> = (0..odd.len()).filter(|&k| idle[k].is_some()).collect();
        if with_idle.len() >= 2 {
            let (k1, k2) = (with_idle[0], with_idle[1]);
            let (q1, q2) = (idle[k1].unwrap(), idle[k2].unwrap());
            for (q, c) in [(q1, odd[k1]), (q2, odd[k2])] {
                core_of[q] = None;
                free[c] += 1;
            }
            gates.push((q1, q2));
        } else if let (Some(&k), true) = (with_idle.first(), odd.len() >= 2) {
            let to = *odd.iter().find(|&&c| c != odd[k]).expect("two odd cores");
            core_of[idle[k].unwrap()] = Some(to);
            free[odd[k]] += 1;
            free[to] -= 1;
        } else {
            return false;
        }
    }
    true
}

/// The idle qubit (no gate in slice `t`) of `core` least attracted to it.
fn least_attached_idle(
    slices: &Timeslices,
    t: usize,
    core: usize,
    core_of: &[Option<usize>],
) -> Option<usize> {
    (0..core_of.len())
        .filter(|&q| core_of[q] == Some(core) && slices.partner(t, q).is_none())
        .min_by(|&x, &y| {
            attraction(slices, t, x, core, core_of)
                .total_cmp(&attraction(slices, t, y, core, core_of))
        })
}

#[cfg(test)]
mod tests {
    use super::super::initial_placement;
    use super::*;
    use crate::rng::Rng;

    fn random_pairs(rng: &mut Rng, qubits: usize, gates: usize) -> Vec<(usize, usize)> {
        (0..gates)
            .map(|_| {
                let a = rng.below(qubits);
                (a, (a + 1 + rng.below(qubits - 1)) % qubits)
            })
            .collect()
    }

    #[test]
    fn mappings_are_valid_full_and_with_spare_room() {
        let mut rng = Rng::new(1);
        // Full and with spare room, even and odd core sizes.
        for (cores, per_core, qubits) in [
            (4, 4, 16),
            (4, 5, 17),
            (10, 10, 100),
            (3, 6, 9),
            (5, 4, 20),
            (2, 3, 6),
            (5, 3, 15),
            (4, 5, 20),
        ] {
            let arch = Architecture::new(cores, per_core).unwrap();
            for _ in 0..5 {
                let slices =
                    Timeslices::from_pairs(qubits, &random_pairs(&mut rng, qubits, 6 * qubits))
                        .unwrap()
                        .fitted_to(&arch);
                let initial = initial_placement(&slices, &arch);
                for lookahead in [true, false] {
                    let options = HqaOptions {
                        lookahead,
                        hops: true,
                    };
                    let mapping = hqa(&slices, &arch, &initial, options).unwrap();
                    mapping.validate(&arch, &slices).unwrap();
                }
            }
        }
    }

    #[test]
    fn one_unfeasible_gate_costs_one_teleportation_when_there_is_room() {
        // Two cores of 3 qubits, 4 qubits: q0 q1 in core 0, q2 q3 in core 1.
        let arch = Architecture::new(2, 3).unwrap();
        let slices = Timeslices::from_pairs(4, &[(1, 2)]).unwrap();
        let mapping = hqa(
            &slices,
            &arch,
            &vec![0, 0, 1, 1],
            HqaOptions {
                lookahead: false,
                hops: true,
            },
        )
        .unwrap();
        mapping.validate(&arch, &slices).unwrap();
        assert_eq!(mapping.communications(), 1);
    }

    #[test]
    fn full_cores_need_a_swap() {
        // Two full cores of 2: a gate across them moves two qubits.
        let arch = Architecture::new(2, 2).unwrap();
        let slices = Timeslices::from_pairs(4, &[(1, 2)]).unwrap();
        let mapping = hqa(
            &slices,
            &arch,
            &vec![0, 0, 1, 1],
            HqaOptions {
                lookahead: false,
                hops: true,
            },
        )
        .unwrap();
        mapping.validate(&arch, &slices).unwrap();
        assert_eq!(mapping.communications(), 2);
    }

    #[test]
    fn attraction_prefers_where_qubits_go_next() {
        // q1–q2 now, then q2–q3: put the pair in q3's core (core 1).
        let arch = Architecture::new(2, 3).unwrap();
        let slices = Timeslices::from_pairs(4, &[(1, 2), (2, 3)]).unwrap();
        let mapping = hqa(&slices, &arch, &vec![0, 0, 1, 1], HqaOptions::default()).unwrap();
        mapping.validate(&arch, &slices).unwrap();
        assert_eq!(mapping.slices[0][1], 1);
        assert_eq!(mapping.communications(), 1);
    }
}
