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
use crate::lower::Lowering;
use crate::passes::{Relabel, best_pipeline, merge_single_qubit_runs, network_promising};
use crate::routing::{CouplingGraph, SabreOptions, layout_and_route};
use crate::two_qubit::kak;
use std::collections::HashSet;

/// SABRE as device compilation uses it: a layout that needs no SWAP when
/// there is one, more layout trials than the routing defaults, and SWAPs
/// that merge with the block before them counted for what they cost (see
/// `routing/sabre.rs`). On the 19-qubit heavy-hex
/// benchmark these took 8% of the two-qubit gates off.
pub fn routing_options() -> SabreOptions {
    SabreOptions {
        trials: 32,
        merge_weight: 0.5,
        perfect_layout: true,
        ..SabreOptions::default()
    }
}

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
    let native: Vec<String> = native.into_iter().collect();
    let graph = CouplingGraph::new(device.qubits, device.edges)?;
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
        let out = on_physical_qubits(&instructions, &native, device)?;
        return Ok(DeviceCircuit {
            instructions: out,
            initial_layout: identity.clone(),
            final_layout: identity,
        });
    }
    // On logical qubits, where relabelling is free (the final layout reports
    // it). A Pauli network can need far more routing than it saves, so a
    // result with one is routed next to one without, and the cheaper kept.
    let mut logical_candidates = vec![best_pipeline(
        &instructions,
        &native,
        Relabel::Always,
        false,
        false,
    )?];
    if network_promising(&instructions, &native) {
        let with_network = best_pipeline(&instructions, &native, Relabel::Always, true, false)?;
        if with_network.0 != logical_candidates[0].0 {
            logical_candidates.push(with_network);
        }
    }
    let mut best: Option<DeviceCircuit> = None;
    for (optimized, moved) in logical_candidates {
        let routed = layout_and_route(&optimized, logical, &graph, options)?;
        let final_layout = (0..logical)
            .map(|q| routed.final_layout[moved.get(q).copied().unwrap_or(q)])
            .collect();
        let out = on_physical_qubits(&routed.instructions, &native, device)?;
        if best
            .as_ref()
            .is_none_or(|b| cost(&out) < cost(&b.instructions))
        {
            best = Some(DeviceCircuit {
                instructions: out,
                initial_layout: routed.initial_layout,
                final_layout,
            });
        }
    }
    Ok(best.expect("a candidate"))
}

/// The pipeline on physical qubits: no relabelling, every pass keeps the
/// pairs; then gates turned around where the couplings have a direction.
fn on_physical_qubits(
    physical: &[Instruction],
    native: &[String],
    device: &Device,
) -> Result<Vec<Instruction>, String> {
    let (mut out, _) = best_pipeline(physical, native, Relabel::Never, false, true)?;
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
                out.push(Instruction::gate(
                    &inst.name,
                    vec![a, b],
                    inst.params.clone(),
                ));
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
                    "cx" | "cz" | "swap" => Instruction::gate(name, vec![a, b], vec![]),
                    "rzz" => Instruction::gate(name, vec![a, b], vec![rng.range(-3.0, 3.0)]),
                    "rz" => Instruction::gate(name, vec![a], vec![rng.range(-3.0, 3.0)]),
                    _ => Instruction::gate(name, vec![a], vec![]),
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
                circuit.push(Instruction::gate("swap", vec![w, x], vec![]));
                let y = at.iter().position(|&v| v == x).expect("a permutation");
                at[y] = w;
                at[x] = x;
            }
        }
        circuit
    }

    fn check(directed: bool, native: &[&str]) {
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
            let out = compile_for_device(
                circuit.clone(),
                native.clone(),
                &device,
                &SabreOptions::default(),
            )
            .unwrap();
            for inst in out.instructions.iter().filter(|i| i.qubits.len() == 2) {
                let pair = (inst.qubits[0], inst.qubits[1]);
                assert!(
                    edges.contains(&pair) || (!directed && edges.contains(&(pair.1, pair.0))),
                    "{inst:?} is not on an allowed pair"
                );
                assert!(native.contains(&inst.name));
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
            Instruction::gate("cx", vec![0, 1], vec![]),
            Instruction::gate("cx", vec![2, 1], vec![]),
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
