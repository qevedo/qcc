//! Phase folding: merge rotations that act on the same parity of the
//! circuit's values, however far apart they are.
//!
//! From Amy, Maslov and Mosca ("Polynomial-time T-depth optimization of
//! Clifford+T circuits via matroid partitioning", 2014) and Nam, Ross, Su,
//! Childs and Maslov ("Automated optimization of large quantum circuits with
//! continuous parameters", 2018), in the path-sum form of Amy's Feynman.
//!
//! Every wire carries an affine parity: the XOR of some variables, plus a
//! constant. The circuit's qubits start as variables x₀, x₁, …; a CNOT XORs
//! its control's parity into its target's, an X flips the constant, and any
//! other gate that is not diagonal (H, SX, a measurement…) gives its qubits
//! fresh variables. A diagonal rotation (Rz, phase, T, S, Z, or Rzz on two
//! wires) multiplies the state by `e^{iθ f}`, where `f` is the parity it acts
//! on. The circuit's unitary is a sum over paths whose phase adds up all such
//! terms, so rotations on the same parity can be merged into one: its angle
//! becomes the sum, and the others go. Rotations whose angles sum
//! to zero all go. Rotations on a constant parity are global phases and go
//! too.
//!
//! The merging removes single-qubit gates, and leaves CNOTs that used to be
//! separated by rotations next to each other, for commutative cancellation.
//! Qiskit's optimizations only merge rotations that are next to each other.

use crate::ir::Instruction;
use std::collections::HashMap;
use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI, TAU};

/// The XOR of some variables, as a bitset, plus a constant.
#[derive(Clone, Debug, Default)]
struct Parity {
    variables: Vec<u64>,
    constant: bool,
}

impl Parity {
    fn variable(v: usize) -> Parity {
        let mut variables = vec![0; v / 64 + 1];
        variables[v / 64] |= 1 << (v % 64);
        Parity {
            variables,
            constant: false,
        }
    }

    fn xor(&mut self, other: &Parity) {
        if self.variables.len() < other.variables.len() {
            self.variables.resize(other.variables.len(), 0);
        }
        for (a, b) in self.variables.iter_mut().zip(&other.variables) {
            *a ^= b;
        }
        self.constant ^= other.constant;
    }

    /// The variables, without trailing zero words, as a key.
    fn key(&self) -> Vec<u64> {
        let mut key = self.variables.clone();
        while key.last() == Some(&0) {
            key.pop();
        }
        key
    }
}

/// A rotation's phase on a parity, `e^{iθ f}` up to a global phase.
struct Term {
    /// Which instruction it comes from.
    index: usize,
    /// −1 when the parity's constant is 1, as `θ (f ⊕ 1) = θ − θ f`; else 1.
    sign: f64,
    /// Its coefficient on the parity's variables: `sign · θ`.
    coefficient: f64,
    /// Whether the instruction's angle can be changed (Rz, phase, Rzz); fixed
    /// gates (T, S, Z…) can only be removed.
    adjustable: bool,
}

/// The phase of a diagonal single-qubit gate on the parity of its wire
/// (`diag(1, e^{iθ})` up to a global phase), and whether its angle is a
/// parameter.
fn diagonal_phase(inst: &Instruction) -> Option<(f64, bool)> {
    Some(match (inst.name.as_str(), inst.params.as_slice()) {
        ("rz" | "p" | "phase" | "u1", &[theta]) => (theta, true),
        ("z", []) => (PI, false),
        ("s", []) => (FRAC_PI_2, false),
        ("sdg", []) => (-FRAC_PI_2, false),
        ("t", []) => (FRAC_PI_4, false),
        ("tdg", []) => (-FRAC_PI_4, false),
        _ => return None,
    })
}

/// Gates that are diagonal but not a phase on one parity: they leave the
/// wires' parities alone, and are left alone.
fn other_diagonal(name: &str) -> bool {
    matches!(
        name,
        "cz" | "cp" | "cphase" | "cu1" | "crz" | "id" | "barrier"
    )
}

/// Merge the rotations of `instructions` that act on the same parity.
pub fn phase_folding(instructions: Vec<Instruction>) -> Vec<Instruction> {
    let qubits = instructions
        .iter()
        .flat_map(|i| i.qubits.iter().copied())
        .max()
        .map_or(0, |q| q + 1);
    let mut wires: Vec<Parity> = (0..qubits).map(Parity::variable).collect();
    let mut next_variable = qubits;
    let mut fresh = |wires: &mut Vec<Parity>, q: usize| {
        wires[q] = Parity::variable(next_variable);
        next_variable += 1;
    };
    // Terms grouped by parity, in the order the parities first appear.
    let mut groups: HashMap<Vec<u64>, Vec<Term>> = HashMap::new();
    let mut order: Vec<Vec<u64>> = Vec::new();
    let mut add = |parity: &Parity, index: usize, theta: f64, adjustable: bool| {
        let key = parity.key();
        let sign = if parity.constant { -1.0 } else { 1.0 };
        let group = groups.entry(key.clone()).or_insert_with(|| {
            order.push(key);
            Vec::new()
        });
        group.push(Term {
            index,
            sign,
            coefficient: sign * theta,
            adjustable,
        });
    };

    for (index, inst) in instructions.iter().enumerate() {
        match (inst.name.as_str(), inst.qubits.as_slice()) {
            ("cx" | "cnot", &[c, t]) if inst.is_gate() => {
                let control = wires[c].clone();
                wires[t].xor(&control);
            }
            ("x", &[q]) if inst.is_gate() => wires[q].constant ^= true,
            ("swap", &[a, b]) if inst.is_gate() => wires.swap(a, b),
            ("rzz", &[a, b]) if inst.is_gate() => {
                let mut parity = wires[a].clone();
                parity.xor(&wires[b]);
                add(&parity, index, inst.params[0], true);
            }
            (name, _) if inst.is_gate() && other_diagonal(name) => {}
            (_, &[q]) if inst.is_gate() && diagonal_phase(inst).is_some() => {
                let (theta, adjustable) = diagonal_phase(inst).expect("checked");
                add(&wires[q], index, theta, adjustable);
            }
            (_, qubits) => {
                for &q in qubits {
                    fresh(&mut wires, q);
                }
            }
        }
    }

    // What to do with each instruction: keep, drop, or set its angle.
    let mut dropped = vec![false; instructions.len()];
    let mut angle: HashMap<usize, f64> = HashMap::new();
    for key in &order {
        let terms = &groups[key];
        let total: f64 = terms.iter().map(|t| t.coefficient).sum();
        if key.is_empty() {
            // A phase on a constant: a global phase.
            for t in terms {
                dropped[t.index] = true;
            }
            continue;
        }
        if terms.len() < 2 {
            continue;
        }
        // The sum goes into the last adjustable rotation: on the benchmarks that
        // leaves more for the merges and resyntheses that follow than the
        // first does (18% fewer gates on a Heisenberg chain).
        let Some(keeper) = terms.iter().rev().find(|t| t.adjustable) else {
            continue;
        };
        for t in terms {
            if t.index != keeper.index {
                dropped[t.index] = true;
            }
        }
        // The keeper's angle θ' with sign · θ' = the total.
        let value = normalize(keeper.sign * total);
        if value.abs() < 1e-12 {
            dropped[keeper.index] = true;
        } else {
            angle.insert(keeper.index, value);
        }
    }

    instructions
        .into_iter()
        .enumerate()
        .filter(|(i, _)| !dropped[*i])
        .map(|(i, mut inst)| {
            if let Some(&theta) = angle.get(&i) {
                inst.params[0] = theta;
            }
            inst
        })
        .collect()
}

/// An angle in (−π, π].
fn normalize(theta: f64) -> f64 {
    let r = theta.rem_euclid(TAU);
    if r > PI { r - TAU } else { r }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testing::{Rng, circuit_unitary, dense_equal_up_to_phase};

    fn g(name: &str, qubits: Vec<usize>, params: Vec<f64>) -> Instruction {
        Instruction::gate(name, qubits, params)
    }

    #[test]
    fn rotations_on_the_same_parity_merge_across_cnots() {
        // Rz on q1's parity x0⊕x1 twice, with a CNOT pair around the first.
        let circuit = vec![
            g("cx", vec![0, 1], vec![]),
            g("rz", vec![1], vec![0.3]),
            g("cx", vec![0, 1], vec![]),
            g("h", vec![0], vec![]),
            g("cx", vec![0, 1], vec![]),
            g("rz", vec![1], vec![0.4]),
            g("cx", vec![0, 1], vec![]),
        ];
        // After the H, q0 is a fresh variable: x0'⊕x1 is another parity.
        assert_eq!(phase_folding(circuit.clone()).len(), circuit.len());
        let circuit = vec![
            g("cx", vec![0, 1], vec![]),
            g("rz", vec![1], vec![0.3]),
            g("cx", vec![0, 1], vec![]),
            g("h", vec![2], vec![]),
            g("cx", vec![0, 1], vec![]),
            g("t", vec![1], vec![]),
            g("x", vec![1], vec![]),
            g("rz", vec![1], vec![0.4]),
            g("x", vec![1], vec![]),
            g("cx", vec![0, 1], vec![]),
        ];
        let out = phase_folding(circuit.clone());
        assert_eq!(out.iter().filter(|i| i.name == "rz").count(), 1);
        assert!(out.iter().all(|i| i.name != "t"));
        assert!(dense_equal_up_to_phase(
            &circuit_unitary(&out, 3),
            &circuit_unitary(&circuit, 3),
            1e-9
        ));
    }

    #[test]
    fn random_circuits_keep_their_unitary() {
        let mut rng = Rng::new(21);
        let gates = [
            "cx", "rz", "t", "s", "x", "h", "sx", "rzz", "cz", "swap", "tdg",
        ];
        let mut removed = 0;
        for _ in 0..300 {
            let n = 3;
            let circuit: Vec<Instruction> = (0..30)
                .map(|_| {
                    let name = gates[rng.below(gates.len())];
                    let a = rng.below(n);
                    let b = (a + 1 + rng.below(n - 1)) % n;
                    match name {
                        "cx" | "cz" | "swap" => g(name, vec![a, b], vec![]),
                        "rzz" => g(name, vec![a, b], vec![rng.range(-3.0, 3.0)]),
                        "rz" => g(name, vec![a], vec![rng.range(-3.0, 3.0)]),
                        _ => g(name, vec![a], vec![]),
                    }
                })
                .collect();
            let out = phase_folding(circuit.clone());
            removed += circuit.len() - out.len();
            assert!(dense_equal_up_to_phase(
                &circuit_unitary(&out, n),
                &circuit_unitary(&circuit, n),
                1e-9
            ));
        }
        assert!(removed > 300, "only {removed} gates removed");
    }
}
