//! SWAP elision: remove SWAP gates by relabelling the qubits of everything
//! that follows them, as Qiskit's `ElidePermutations` does.
//!
//! A SWAP costs three CNOTs; relabelling costs nothing, but leaves the
//! qubits' states on other wires at the end. That cannot be observed when
//! each moved state is measured last, and is otherwise the caller's to
//! accept: `wires` says where each state ends.

use crate::ir::Instruction;

/// Remove the SWAP gates of `instructions`, on `qubits` qubits, relabelling
/// what follows each. Returns the instructions and `wires`: qubit `q`'s state
/// ends on wire `wires[q]`.
pub fn elide_swaps(
    instructions: Vec<Instruction>,
    qubits: usize,
) -> (Vec<Instruction>, Vec<usize>) {
    let qubits = instructions
        .iter()
        .flat_map(|i| i.qubits.iter().map(|&q| q + 1))
        .fold(qubits, usize::max);
    // The wire that holds each qubit's state.
    let mut wires: Vec<usize> = (0..qubits).collect();
    let mut out = Vec::with_capacity(instructions.len());
    for inst in instructions {
        match inst.qubits.as_slice() {
            &[a, b] if inst.name == "swap" && inst.is_gate() => wires.swap(a, b),
            _ => out.push(inst.relabel(&wires)),
        }
    }
    (out, wires)
}

/// Whether the relabelling `wires` cannot be observed: every wire that ends
/// up holding another qubit's state is measured last (barriers aside).
pub fn measured_out(instructions: &[Instruction], wires: &[usize]) -> bool {
    let mut last: Vec<Option<&str>> = vec![None; wires.len()];
    for inst in instructions.iter().filter(|i| i.name != "barrier") {
        for &q in &inst.qubits {
            last[q] = Some(inst.name.as_str());
        }
    }
    wires
        .iter()
        .enumerate()
        .all(|(q, &w)| w == q || last[w] == Some("measure"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testing::{Rng, circuit_unitary, dense_equal_up_to_phase};

    fn g(name: &str, qubits: Vec<usize>) -> Instruction {
        Instruction::gate(name, qubits, vec![])
    }

    fn measure(q: usize) -> Instruction {
        Instruction {
            name: "measure".into(),
            qubits: vec![q],
            params: vec![],
            clbits: vec![q],
        }
    }

    #[test]
    fn swaps_become_a_relabelling() {
        let circuit = vec![g("h", vec![0]), g("swap", vec![0, 1]), g("cx", vec![1, 2])];
        let (out, wires) = elide_swaps(circuit, 3);
        assert_eq!(out, vec![g("h", vec![0]), g("cx", vec![0, 2])]);
        assert_eq!(wires, vec![1, 0, 2]);
        assert!(!measured_out(&out, &wires));
        let mut measured = out.clone();
        measured.extend([measure(0), measure(1)]);
        assert!(measured_out(&measured, &wires));
    }

    #[test]
    fn random_circuits_agree_up_to_the_permutation() {
        let mut rng = Rng::new(5);
        let names = ["h", "t", "cx", "swap", "swap"];
        for _ in 0..100 {
            let n = 3;
            let circuit: Vec<Instruction> = (0..20)
                .map(|_| {
                    let name = names[rng.below(names.len())];
                    let a = rng.below(n);
                    let b = (a + 1 + rng.below(n - 1)) % n;
                    match name {
                        "cx" | "swap" => g(name, vec![a, b]),
                        _ => g(name, vec![a]),
                    }
                })
                .collect();
            let (mut out, mut wires) = elide_swaps(circuit.clone(), n);
            assert!(out.iter().all(|i| i.name != "swap"));
            // Moving each state back to its own wire restores the unitary.
            for q in 0..n {
                let w = wires[q];
                if w != q {
                    out.push(g("swap", vec![w, q]));
                    let p = wires.iter().position(|&x| x == q).expect("a permutation");
                    wires[p] = w;
                    wires[q] = q;
                }
            }
            assert!(dense_equal_up_to_phase(
                &circuit_unitary(&out, n),
                &circuit_unitary(&circuit, n),
                1e-9
            ));
        }
    }
}
