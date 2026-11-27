//! Compiling for a device whose qubits are not all coupled.
//!
//! The circuit is optimized on logical qubits first (fewer two-qubit gates
//! leave the router less to do), laid out and routed with SABRE
//! ([`crate::routing`]), and optimized again on the physical qubits, where
//! every pass keeps each two-qubit gate on the pair it was routed to. On a
//! device whose couplings have a direction, gates the other way round are
//! turned around with single-qubit gates. A circuit whose two-qubit gates
//! all act on coupled pairs already is not routed again.

use crate::gates::matrix_2q;
use crate::ir::{Instruction, cost};
use crate::linalg::{M2, equal_up_to_phase4, kron};
use crate::lower::{Lowering, SynthesisCache};
use crate::passes::{
    Last, commutative_cancellation, lowered_for_network, merge_single_qubit_runs, network_trials,
    on_the_same_pairs, with_and_without_network,
};
use crate::pauli_network::{merge_rotations, placement, synthesize_on};
use crate::phase_folding::phase_folding;
use crate::routing::{CouplingGraph, SabreOptions, layout_and_route_and_first, route};
use crate::two_qubit::kak;
use std::collections::HashSet;

/// SABRE as device compilation uses it: a layout that needs no SWAP when
/// there is one, more layout trials than the routing defaults, and SWAPs
/// that merge with the block before them counted for what they cost (see
/// `routing/sabre.rs`). On the 19-qubit heavy-hex
/// benchmark these took 8% of the two-qubit gates off.
///
/// Each layout is also routed 16 times with different tie-breaking rather
/// than the default 4 ([`USUAL_SWAP_TRIALS`]): SABRE's choices vary a lot
/// from one tie-breaking to the next (Qiskit's own level 3 leaves Quantum
/// Volume at 16 qubits with 798 to 840 CNOTs on this chip, depending on its
/// seed). Routing itself got 2–3× cheaper in the same change (see
/// `routing/sabre.rs`), so the wider search costs about what four routings
/// per layout used to.
pub fn routing_options() -> SabreOptions {
    SabreOptions {
        trials: 32,
        swap_trials: 16,
        merge_weight: 0.5,
        perfect_layout: true,
        ..SabreOptions::default()
    }
}

/// The routings per layout of SABRE's defaults, and of device compilation
/// before it searched more widely: the routing those would have chosen is
/// optimized too, next to the wider search's (see [`compile_for_device`]).
pub const USUAL_SWAP_TRIALS: usize = 4;

/// The device: its qubits and the pairs it couples.
pub struct Device<'a> {
    pub qubits: usize,
    pub edges: &'a [(usize, usize)],
    /// Whether a two-qubit gate runs only from the first qubit of an edge to
    /// its second.
    pub directed: bool,
}

/// A circuit compiled for a device.
#[derive(Clone, Debug)]
pub struct DeviceCircuit {
    /// The circuit, on the device's physical qubits.
    pub instructions: Vec<Instruction>,
    /// `initial_layout[q]`: the physical qubit that holds logical qubit `q` at
    /// the start.
    pub initial_layout: Vec<usize>,
    /// `final_layout[q]`: the physical qubit that holds its state at the end.
    pub final_layout: Vec<usize>,
}

/// Compile `instructions` for `device` and its `native` gates.
pub fn compile_for_device(
    instructions: Vec<Instruction>,
    native: impl IntoIterator<Item = String>,
    device: &Device,
    options: &SabreOptions,
) -> Result<DeviceCircuit, String> {
    compile_with(instructions, native, device, options, true)
}

/// [`compile_for_device`], with or without the Pauli networks synthesized
/// on the chip (`on_chip`; the tests compare the two).
fn compile_with(
    instructions: Vec<Instruction>,
    native: impl IntoIterator<Item = String>,
    device: &Device,
    options: &SabreOptions,
    on_chip: bool,
) -> Result<DeviceCircuit, String> {
    let native: Vec<String> = native.into_iter().collect();
    let graph = CouplingGraph::new(device.qubits, device.edges)?;
    // One cache of two-qubit syntheses for every pipeline below: on logical
    // qubits, then on each candidate's physical qubits, round after round.
    let cache = SynthesisCache::new();
    let logical = instructions
        .iter()
        .flat_map(|i| i.qubits.iter().map(|&q| q + 1))
        .max()
        .unwrap_or(0);
    if logical > device.qubits {
        return Err(format!(
            "{logical} qubits do not fit on a device of {}",
            device.qubits
        ));
    }
    let coupled = |a: usize, b: usize| graph.distance(a, b) == 1;
    let routed_already = instructions
        .iter()
        .filter(|i| i.is_gate() && i.qubits.len() == 2)
        .all(|i| coupled(i.qubits[0], i.qubits[1]));
    if routed_already {
        let identity: Vec<usize> = (0..logical).collect();
        let out = on_physical_qubits(&instructions, &native, device, &cache)?;
        return Ok(DeviceCircuit {
            instructions: out,
            initial_layout: identity.clone(),
            final_layout: identity,
        });
    }
    // On logical qubits, where relabelling is free (the final layout reports
    // it). A Pauli network can need far more routing than it saves, so a
    // result with one is routed next to one without, and the cheaper kept.
    let logical_candidates = with_and_without_network(&instructions, &native, Last::No, &cache)?;
    // Each candidate is laid out and routed on its own thread. The CNOTs a
    // routing's SWAPs cost are only an estimate of what the optimized circuit
    // will need (the single-qubit gates, a CZ device's rounds and a Pauli
    // network's local Cliffords all move it), so the cheapest routing of the
    // wide search can come out worse than the one the usual four routings
    // per layout would have kept: both are optimized, and the cheaper kept.
    // On 40 held-out circuits on the benchmark chip, the wide search alone
    // left 3 with more two-qubit gates than before; with both, none.
    let routings: Vec<(usize, crate::routing::Routed)> =
        crate::parallel::map_ordered(logical_candidates.len(), |i| {
            layout_and_route_and_first(
                &logical_candidates[i].0,
                logical,
                &graph,
                options,
                USUAL_SWAP_TRIALS,
            )
        })
        .into_iter()
        .enumerate()
        .map(|(i, routed)| routed.map(|r| r.into_iter().map(move |r| (i, r))))
        .collect::<Result<Vec<_>, String>>()?
        .into_iter()
        .flatten()
        .collect();
    // What each routing leaves: the routed circuit, its initial layout and
    // where each logical qubit's state ends.
    let mut jobs: Vec<DeviceCircuit> = routings
        .iter()
        .map(|(i, routed)| {
            let (_, moved) = &logical_candidates[*i];
            let final_layout: Vec<usize> = (0..logical)
                .map(|q| routed.final_layout[moved.get(q).copied().unwrap_or(q)])
                .collect();
            DeviceCircuit {
                instructions: routed.instructions.clone(),
                initial_layout: routed.initial_layout.clone(),
                final_layout,
            }
        })
        .collect();
    // A circuit of Pauli rotations is also synthesized again on the chip's
    // own pairs, from the layouts the routings found, and routed as one more
    // candidate (see `networks_on_chip`).
    let lowered = on_chip
        .then(|| lowered_for_network(&instructions, &native))
        .flatten();
    if let Some(lowered) = lowered {
        jobs.extend(networks_on_chip(&lowered, logical, &jobs, &graph, options));
    }
    // Each routing is then optimized on the chip's qubits on its own thread
    // (the rounds of `until_fixed_point` make that most of the time for CZ),
    // and the cheapest kept: the first on ties, as in a loop.
    let routed = crate::parallel::map_ordered(jobs.len(), |j| {
        let job = &jobs[j];
        let out = on_physical_qubits(&job.instructions, &native, device, &cache)?;
        Ok::<_, String>(DeviceCircuit {
            instructions: out,
            initial_layout: job.initial_layout.clone(),
            final_layout: job.final_layout.clone(),
        })
    });
    let mut best: Option<DeviceCircuit> = None;
    for candidate in routed {
        let candidate = candidate?;
        if best
            .as_ref()
            .is_none_or(|b| cost(&candidate.instructions) < cost(&b.instructions))
        {
            best = Some(candidate);
        }
    }
    Ok(best.expect("a candidate"))
}

/// A Pauli network synthesized on the chip's pairs, as a routed candidate:
/// for each initial layout of `routed` (the routings of the other
/// candidates), each placement of the rotations improved from one of them
/// or grown greedily (`placement` in `pauli_network.rs`), and each of the
/// pipeline's network trials, the rotations of
/// `lowered` placed by that layout and synthesized again with every CX that
/// sets them up on a coupled pair (`synthesize_on` in `pauli_network.rs`),
/// each with two endings: the Clifford left at the end synthesized on the
/// chip's pairs too, or synthesized as off a device, whose CXs are not, and
/// routed from there. Of each ending, the trial needing fewest CXs and
/// SWAPs once its cancelling CXs and mergeable rotations are gone, which is
/// then optimized like the others. The Clifford on the chip wins the
/// ranking on nearly every trial, but not every optimized result: with it
/// alone (in a first version, which gathered each pair's X image first),
/// the 8-qubit benchmark Hamiltonian got 191 CNOTs instead of 0.4.15's 189.
///
/// Synthesized off the chip and routed, a network's CXs land on pairs far
/// apart: the 16-qubit benchmark Hamiltonian needs 313 CXs on all pairs and
/// 1321 routed on the heavy-hex chip. The greedy synthesis on the chip
/// leaves many CXs that cancel (about 3000 raw for the same Hamiltonian),
/// so the candidates are ranked after a light cleanup: ranked on raw
/// counts, the one optimized was rarely the best.
fn networks_on_chip(
    lowered: &[Instruction],
    logical: usize,
    routed: &[DeviceCircuit],
    graph: &CouplingGraph,
    options: &SabreOptions,
) -> Vec<DeviceCircuit> {
    let mut layouts: Vec<Vec<usize>> = Vec::new();
    for r in routed {
        if !layouts.contains(&r.initial_layout) {
            layouts.push(r.initial_layout.clone());
        }
    }
    // Besides SABRE's layouts, placements chosen for the rotations
    // themselves (`placement`): one grown greedily, and each of SABRE's
    // improved by exchanges. SABRE lays out the circuit as written, CX
    // ladders and all, and its layouts are not the best for the network:
    // on 22 held-out Hamiltonians (6–16 qubits) on the benchmark chip, these
    // took the CNOTs from 10463 to 10068, more on none. Without SABRE's own
    // layouts they would leave 10265 (8 and 16 qubits got worse), and with
    // only half the trials on the new ones 10142.
    let from_sabre = layouts.len();
    let mut placed = Vec::new();
    placed.extend(placement(lowered, logical, graph, None));
    for layout in &layouts {
        placed.extend(placement(lowered, logical, graph, Some(layout)));
    }
    for layout in placed {
        if !layouts.contains(&layout) {
            layouts.push(layout);
        }
    }
    let n = graph.qubits();
    let identity: Vec<usize> = (0..n).collect();
    let trials: Vec<(usize, &Vec<usize>, (u64, f64))> = layouts
        .iter()
        .enumerate()
        .flat_map(|(l, layout)| {
            let placed = usize::from(l >= from_sabre);
            network_trials().map(move |trial| (placed, layout, trial))
        })
        .collect();
    let candidates = crate::parallel::map_ordered(trials.len(), |t| {
        let (_, layout, (seed, penalty)) = trials[t];
        let placed: Vec<Instruction> = lowered.iter().map(|i| i.relabel(layout)).collect();
        let networks = synthesize_on(&placed, graph, seed, penalty)?;
        let mut kinds = Vec::new();
        for (network, wires) in networks {
            let network = commutative_cancellation(phase_folding(merge_rotations(network)));
            let Ok(routed) = route(&network, n, graph, &identity, options) else {
                kinds.push(None);
                continue;
            };
            let cxs = routed
                .instructions
                .iter()
                .filter(|i| i.name == "cx")
                .count();
            // Logical qubit q starts on layout[q], the network leaves that
            // state on wires[layout[q]], and routing moves it from there.
            let final_layout = (0..logical)
                .map(|q| routed.final_layout[wires[layout[q]]])
                .collect();
            let circuit = DeviceCircuit {
                instructions: routed.instructions,
                initial_layout: layout.clone(),
                final_layout,
            };
            kinds.push(Some((cxs + routed.swap_cost, circuit)));
        }
        Some(kinds)
    });
    // Of each kind of tail, the first of the cheapest, as in a loop; from
    // SABRE's layouts and from the placements apart, so that the candidates
    // of SABRE's layouts are still optimized: a placement can win the
    // ranking and lose once optimized (with one best of all, a held-out
    // 9-qubit Hamiltonian with CZ got 259 CZs instead of 258).
    let mut best: [Vec<Option<(usize, DeviceCircuit)>>; 2] = [Vec::new(), Vec::new()];
    for (&(placed, _, _), kinds) in trials.iter().zip(candidates) {
        let Some(kinds) = kinds else { continue };
        let best = &mut best[placed];
        best.resize_with(best.len().max(kinds.len()), || None);
        for (b, candidate) in best.iter_mut().zip(kinds) {
            let Some(candidate) = candidate else { continue };
            if b.as_ref().is_none_or(|(c, _)| candidate.0 < *c) {
                *b = Some(candidate);
            }
        }
    }
    best.into_iter()
        .flatten()
        .flatten()
        .map(|(_, circuit)| circuit)
        .collect()
}

/// The pipeline on physical qubits: no relabelling, every pass keeps the
/// pairs, repeated while it shrinks the circuit (for CZ-like devices); then
/// gates turned around where the couplings have a direction.
fn on_physical_qubits(
    physical: &[Instruction],
    native: &[String],
    device: &Device,
    cache: &SynthesisCache,
) -> Result<Vec<Instruction>, String> {
    let mut out = on_the_same_pairs(physical, native, cache)?;
    if device.directed {
        let allowed: HashSet<(usize, usize)> = device.edges.iter().copied().collect();
        out = turn_around(out, &allowed, native);
    }
    Ok(out)
}

/// Rewrite two-qubit gates on `(b, a)` where only `(a, b)` is allowed as
/// single-qubit gates around the gate on `(a, b)`, found from the KAK
/// decompositions of the gate both ways round, then merge the runs.
fn turn_around(
    instructions: Vec<Instruction>,
    allowed: &HashSet<(usize, usize)>,
    native: &[String],
) -> Vec<Instruction> {
    let lowering = Lowering::new(native.iter().cloned());
    let mut out = Vec::with_capacity(instructions.len());
    for inst in instructions {
        let reversible = inst.is_gate() && inst.qubits.len() == 2;
        let (b, a) = match inst.qubits.as_slice() {
            &[b, a] if reversible && !allowed.contains(&(b, a)) && allowed.contains(&(a, b)) => {
                (b, a)
            }
            _ => {
                out.push(inst);
                continue;
            }
        };
        match reversal(&inst.name, &inst.params) {
            Some([before_a, before_b, after_a, after_b]) => {
                for (q, m) in [(a, before_a), (b, before_b)] {
                    out.extend(lowering.one_qubit(&m, q).unwrap_or_default());
                }
                out.push(Instruction::gate(&inst.name, [a, b], inst.params.clone()));
                for (q, m) in [(a, after_a), (b, after_b)] {
                    out.extend(lowering.one_qubit(&m, q).unwrap_or_default());
                }
            }
            None => out.push(inst),
        }
    }
    merge_single_qubit_runs(out, &lowering.native, &lowering.bases_1q)
}

/// Single-qubit gates `[before_a, before_b, after_a, after_b]` with
/// g on (b, a) = (after_a ⊗ after_b) · g on (a, b) · (before_a ⊗ before_b),
/// `a` the first factor; None if the gate is unknown.
fn reversal(name: &str, params: &[f64]) -> Option<[M2; 4]> {
    let forward = matrix_2q(name, params)?;
    let swap = matrix_2q("swap", &[])?;
    let reversed = swap * forward * swap;
    let (f, r) = (kak(&forward), kak(&reversed));
    // forward ∝ (fa1 ⊗ fa0) C (fb1 ⊗ fb0), reversed ∝ (ra1 ⊗ ra0) C (rb1 ⊗ rb0).
    let out = [
        f.b1.adjoint() * r.b1,
        f.b0.adjoint() * r.b0,
        r.a1 * f.a1.adjoint(),
        r.a0 * f.a0.adjoint(),
    ];
    let rebuilt = kron(&out[2], &out[3]) * forward * kron(&out[0], &out[1]);
    equal_up_to_phase4(&rebuilt, &reversed, 1e-9).then_some(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testing::{Rng, circuit_unitary, dense_equal_up_to_phase};

    const INPUT: [&str; 8] = ["h", "t", "rz", "sx", "cx", "cz", "swap", "rzz"];

    fn random_circuit(rng: &mut Rng, n: usize, length: usize) -> Vec<Instruction> {
        (0..length)
            .map(|_| {
                let name = INPUT[rng.below(INPUT.len())];
                let a = rng.below(n);
                let b = (a + 1 + rng.below(n - 1)) % n;
                match name {
                    "cx" | "cz" | "swap" => Instruction::gate(name, [a, b], []),
                    "rzz" => Instruction::gate(name, [a, b], [rng.range(-3.0, 3.0)]),
                    "rz" => Instruction::gate(name, [a], [rng.range(-3.0, 3.0)]),
                    _ => Instruction::gate(name, [a], []),
                }
            })
            .collect()
    }

    /// `circuit` followed by SWAPs that move each state from `from[x]` to `x`.
    fn moved_back(mut circuit: Vec<Instruction>, from: &[usize]) -> Vec<Instruction> {
        let mut at = from.to_vec(); // at[x]: where the state for x is now
        for x in 0..at.len() {
            let w = at[x];
            if w != x {
                circuit.push(Instruction::gate("swap", [w, x], []));
                let y = at.iter().position(|&v| v == x).expect("a permutation");
                at[y] = w;
                at[x] = x;
            }
        }
        circuit
    }

    fn check(directed: bool, native: &[&str]) {
        check_with(directed, native, &SabreOptions::default());
    }

    fn check_with(directed: bool, native: &[&str], options: &SabreOptions) {
        let n = 5;
        let edges: Vec<(usize, usize)> = (0..n - 1).map(|i| (i, i + 1)).collect();
        let device = Device {
            qubits: n,
            edges: &edges,
            directed,
        };
        let native: Vec<String> = native.iter().map(|g| g.to_string()).collect();
        let mut rng = Rng::new(31);
        for _ in 0..20 {
            let circuit = random_circuit(&mut rng, n, 30);
            let out =
                compile_for_device(circuit.clone(), native.clone(), &device, options).unwrap();
            for inst in out.instructions.iter().filter(|i| i.qubits.len() == 2) {
                let pair = (inst.qubits[0], inst.qubits[1]);
                assert!(
                    edges.contains(&pair) || (!directed && edges.contains(&(pair.1, pair.0))),
                    "{inst:?} is not on an allowed pair"
                );
                assert!(native.iter().any(|g| inst.name == *g));
            }
            // The input placed by the initial layout equals the output with
            // every state moved back from the final layout to the initial one.
            let placed: Vec<Instruction> = circuit
                .iter()
                .map(|i| i.relabel(&out.initial_layout))
                .collect();
            let mut from = vec![0; n];
            for q in 0..n {
                from[out.initial_layout[q]] = out.final_layout[q];
            }
            assert!(dense_equal_up_to_phase(
                &circuit_unitary(&moved_back(out.instructions.clone(), &from), n),
                &circuit_unitary(&placed, n),
                1e-7
            ));
        }
    }

    #[test]
    fn routed_circuits_agree_up_to_their_layouts() {
        check(false, &["rz", "sx", "x", "cx"]);
        check(false, &["rz", "sx", "x", "cz"]);
        // With the device options: the wide search, and the usual routing
        // optimized beside it.
        let options = SabreOptions {
            trials: 4,
            ..routing_options()
        };
        check_with(false, &["rz", "sx", "x", "cx"], &options);
        check_with(true, &["rz", "sx", "x", "cz"], &options);
    }

    #[test]
    fn the_wide_search_is_never_worse_than_the_usual_one() {
        // The routing the usual swap trials would have kept is optimized as
        // well, so searching more widely can only make a circuit cheaper.
        let edges = [
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (1, 6),
            (6, 7),
            (7, 8),
            (3, 8),
        ];
        let device = Device {
            qubits: 9,
            edges: &edges,
            directed: false,
        };
        let mut rng = Rng::new(35);
        let (mut cheaper, mut dearer) = (0, 0);
        for native in [["rz", "sx", "x", "cx"], ["rz", "sx", "x", "cz"]] {
            let native: Vec<String> = native.iter().map(|g| g.to_string()).collect();
            for seed in 0..6 {
                let circuit = random_circuit(&mut rng, 9, 150);
                let wide = SabreOptions {
                    trials: 4,
                    seed,
                    ..routing_options()
                };
                let usual = SabreOptions {
                    swap_trials: USUAL_SWAP_TRIALS,
                    ..wide
                };
                let compile = |options: &SabreOptions| {
                    let out = compile_for_device(circuit.clone(), native.clone(), &device, options);
                    cost(&out.unwrap().instructions)
                };
                let (w, u) = (compile(&wide), compile(&usual));
                assert!(w <= u, "{w:?} > {u:?}");
                cheaper += usize::from(w < u);
                dearer += usize::from(w > u);
            }
        }
        assert!(cheaper > 0 && dearer == 0, "{cheaper} {dearer}");
    }

    #[test]
    fn networks_on_the_chip_are_never_worse_and_keep_the_circuit() {
        // A heavy-hex cell (a ring of 6) with two tails: degree 3 at most.
        let edges = [
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (5, 0),
            (1, 6),
            (4, 7),
        ];
        let device = Device {
            qubits: 8,
            edges: &edges,
            directed: false,
        };
        let options = SabreOptions {
            trials: 4,
            ..routing_options()
        };
        let mut rng = Rng::new(41);
        let mut cheaper = 0;
        for native in [["rz", "sx", "x", "cx"], ["rz", "sx", "x", "cz"]] {
            let native: Vec<String> = native.iter().map(|g| g.to_string()).collect();
            for k in 0..4 {
                let n = 6 + k % 3;
                let circuit = crate::testing::pauli_exponentials(&mut rng, n, 3 * n);
                let compile = |on_chip| {
                    compile_with(circuit.clone(), native.clone(), &device, &options, on_chip)
                        .unwrap()
                };
                let (with, without) = (compile(true), compile(false));
                let (w, wo) = (cost(&with.instructions), cost(&without.instructions));
                assert!(w <= wo, "{w:?} > {wo:?}");
                cheaper += usize::from(w.0 < wo.0);
                for inst in with.instructions.iter().filter(|i| i.qubits.len() == 2) {
                    let (a, b) = (inst.qubits[0], inst.qubits[1]);
                    assert!(
                        edges.contains(&(a, b)) || edges.contains(&(b, a)),
                        "{inst:?}"
                    );
                }
                // The input placed by the initial layout equals the output
                // with every state moved back to where it started.
                let placed: Vec<Instruction> = circuit
                    .iter()
                    .map(|i| i.relabel(&with.initial_layout))
                    .collect();
                let mut from: Vec<usize> = (0..device.qubits).collect();
                let mut used = vec![false; device.qubits];
                for q in 0..n {
                    from[with.initial_layout[q]] = with.final_layout[q];
                    used[with.final_layout[q]] = true;
                }
                // Qubits the circuit does not use start in |0⟩ and may be
                // left anywhere; give them the free places in order.
                let mut free = (0..device.qubits).filter(|&p| !used[p]);
                for p in 0..device.qubits {
                    if !with.initial_layout[..n].contains(&p) {
                        from[p] = free.next().expect("a free qubit");
                    }
                }
                assert!(dense_equal_up_to_phase(
                    &circuit_unitary(&moved_back(with.instructions.clone(), &from), 8),
                    &circuit_unitary(&placed, 8),
                    1e-7
                ));
            }
        }
        // Pauli exponentials on random qubits are the kind of circuit the
        // chip's own networks are for: most get fewer two-qubit gates.
        assert!(cheaper >= 4, "{cheaper}");
    }

    #[test]
    fn directed_devices_get_gates_their_way_round() {
        check(true, &["rz", "sx", "x", "cx"]);
        check(true, &["rz", "sx", "x", "ecr"]);
    }

    #[test]
    fn routed_circuits_keep_their_layout() {
        let edges = [(0, 1), (1, 2)];
        let device = Device {
            qubits: 3,
            edges: &edges,
            directed: false,
        };
        let circuit = vec![
            Instruction::gate("cx", [0, 1], []),
            Instruction::gate("cx", [2, 1], []),
        ];
        let native = ["rz", "sx", "x", "cx"].map(String::from);
        let out = compile_for_device(circuit, native, &device, &SabreOptions::default()).unwrap();
        assert_eq!(
            (out.initial_layout, out.final_layout),
            (vec![0, 1, 2], vec![0, 1, 2])
        );
        assert_eq!(out.instructions.len(), 2);
    }
}
