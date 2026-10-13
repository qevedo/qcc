//! Decompose one gate or unitary into a native gate set and check the result:
//! the entry point of the decomposition playground on qevedo.com.

use crate::gates::{matrix_1q, matrix_2q, signature};
use crate::ir::Instruction;
use crate::linalg::{C, M2, M4, ONE, ZERO};
use crate::lower::Lowering;
use crate::passes::{commutative_cancellation, compile, merge_single_qubit_runs};

/// A dense matrix, rows of columns, qubit 0 the most significant.
pub type Dense = Vec<Vec<C>>;

/// A decomposition and how far its unitary is from the input's, after
/// removing the global phase (the largest element-wise difference).
#[derive(Debug)]
pub struct Decomposition {
    pub qubits: usize,
    pub instructions: Vec<Instruction>,
    pub error: f64,
}

/// Decompose standard gate `name(params)` into the `native` gates.
pub fn decompose_gate(
    name: &str,
    params: &[f64],
    native: &[String],
) -> Result<Decomposition, String> {
    let name = name.to_lowercase();
    let (qubits, count) = signature(&name).ok_or_else(|| format!("unknown gate '{name}'"))?;
    if params.len() != count {
        let s = if count == 1 { "" } else { "s" };
        return Err(format!(
            "the gate {name} takes {count} parameter{s}, not {}",
            params.len()
        ));
    }
    let target = gate_dense(&name, params).ok_or_else(|| format!("no matrix for '{name}'"))?;
    let gate = Instruction::gate(&name, (0..qubits).collect(), params.to_vec());
    let instructions = compile(vec![gate], native.iter().cloned())?;
    let error = distance(&unitary(&instructions, qubits), &target);
    Ok(Decomposition {
        qubits,
        instructions,
        error,
    })
}

/// Decompose a 2×2 or 4×4 unitary into the `native` gates.
pub fn decompose_unitary(matrix: &Dense, native: &[String]) -> Result<Decomposition, String> {
    let dim = matrix.len();
    if !(dim == 2 || dim == 4) || matrix.iter().any(|row| row.len() != dim) {
        return Err("the matrix must be 2×2 (one qubit) or 4×4 (two qubits)".into());
    }
    if !is_unitary(matrix) {
        return Err("the matrix is not unitary".into());
    }
    let mut lowering = Lowering::new(native.iter().cloned());
    let instructions = if dim == 2 {
        lowering.one_qubit(
            &M2([[matrix[0][0], matrix[0][1]], [matrix[1][0], matrix[1][1]]]),
            0,
        )?
    } else {
        let mut m = M4([[ZERO; 4]; 4]);
        for i in 0..4 {
            for j in 0..4 {
                m.0[i][j] = matrix[i][j];
            }
        }
        let synthesized = lowering
            .two_qubit(&m, [0, 1])
            .map_err(|e| e.replace("cannot lower ''", "cannot lower the unitary"))?;
        let merged = merge_single_qubit_runs(synthesized, &lowering.native, &lowering.bases_1q);
        merge_single_qubit_runs(
            commutative_cancellation(merged),
            &lowering.native,
            &lowering.bases_1q,
        )
    };
    let qubits = if dim == 2 { 1 } else { 2 };
    let error = distance(&unitary(&instructions, qubits), matrix);
    Ok(Decomposition {
        qubits,
        instructions,
        error,
    })
}

/// The matrix of a standard gate on up to three qubits.
pub fn gate_dense(name: &str, params: &[f64]) -> Option<Dense> {
    if let Some(m) = matrix_1q(name, params) {
        return Some(m.0.iter().map(|r| r.to_vec()).collect());
    }
    if let Some(m) = matrix_2q(name, params) {
        return Some(m.0.iter().map(|r| r.to_vec()).collect());
    }
    // Three-qubit permutations: |abc⟩ → |perm(abc)⟩.
    let perm: fn(usize) -> usize = match (name, params) {
        ("ccx", []) => |i| if i >= 6 { i ^ 1 } else { i },
        ("cswap", []) => |i| if i == 5 || i == 6 { 11 - i } else { i },
        _ => return None,
    };
    let mut m = vec![vec![ZERO; 8]; 8];
    for i in 0..8 {
        m[perm(i)][i] = ONE;
    }
    Some(m)
}

/// The unitary of `instructions` on `n` qubits, qubit 0 the most significant.
pub fn unitary(instructions: &[Instruction], n: usize) -> Dense {
    let dim = 1 << n;
    let mut u: Dense = (0..dim)
        .map(|i| (0..dim).map(|j| if i == j { ONE } else { ZERO }).collect())
        .collect();
    let bit = |index: usize, q: usize| (index >> (n - 1 - q)) & 1;
    for inst in instructions {
        let Some(g) = gate_dense(&inst.name, &inst.params) else {
            continue;
        };
        let k = inst.qubits.len();
        let mut next = vec![vec![ZERO; dim]; dim];
        for row in 0..dim {
            let sub_row = inst
                .qubits
                .iter()
                .fold(0, |acc, &q| (acc << 1) | bit(row, q));
            for (sub_col, &coefficient) in g[sub_row].iter().enumerate() {
                if coefficient == ZERO {
                    continue;
                }
                let mut source = row;
                for (pos, &q) in inst.qubits.iter().enumerate() {
                    let mask = 1 << (n - 1 - q);
                    source = if (sub_col >> (k - 1 - pos)) & 1 == 1 {
                        source | mask
                    } else {
                        source & !mask
                    };
                }
                for col in 0..dim {
                    next[row][col] += coefficient * u[source][col];
                }
            }
        }
        u = next;
    }
    u
}

/// The largest element-wise difference between `a` and `b` after aligning their global phases.
pub fn distance(a: &Dense, b: &Dense) -> f64 {
    let overlap: C = a
        .iter()
        .zip(b)
        .flat_map(|(ra, rb)| ra.iter().zip(rb).map(|(x, y)| y.conj() * x))
        .sum();
    if overlap.norm() < 1e-12 {
        return f64::INFINITY;
    }
    let phase = overlap / overlap.norm();
    a.iter()
        .zip(b)
        .flat_map(|(ra, rb)| ra.iter().zip(rb).map(move |(x, y)| (x - y * phase).norm()))
        .fold(0.0, f64::max)
}

fn is_unitary(m: &Dense) -> bool {
    let n = m.len();
    (0..n).all(|i| {
        (0..n).all(|j| {
            let dot: C = (0..n).map(|k| m[k][i].conj() * m[k][j]).sum();
            (dot - if i == j { ONE } else { ZERO }).norm() < 1e-6
        })
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testing::{RandomUnitary, Rng};

    fn native(names: &[&str]) -> Vec<String> {
        names.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn every_standard_gate_decomposes_exactly() {
        let names = [
            "id", "x", "y", "z", "h", "s", "sdg", "t", "tdg", "sx", "sxdg", "rx", "ry", "rz", "p",
            "u1", "u2", "u3", "u", "cx", "cy", "cz", "ch", "csx", "swap", "iswap", "dcx", "ecr",
            "cp", "crx", "cry", "crz", "rxx", "ryy", "rzz", "rzx", "cu3", "cu", "ccx", "cswap",
        ];
        for set in [
            &["rz", "sx", "x", "cx"][..],
            &["rx", "ry", "rxx"],
            &["u3", "ecr"],
            &["rz", "rx", "cz"],
        ] {
            for name in names {
                let count = signature(name).unwrap().1;
                let params: Vec<f64> = (0..count).map(|i| 0.3 + 0.7 * i as f64).collect();
                let d = decompose_gate(name, &params, &native(set)).unwrap();
                assert!(d.error < 1e-9, "{name} on {set:?}: {}", d.error);
                assert!(
                    d.instructions
                        .iter()
                        .all(|i| set.contains(&i.name.as_str()))
                );
            }
        }
    }

    #[test]
    fn unitaries_decompose_exactly() {
        let mut rng = Rng::new(1);
        for _ in 0..50 {
            let u2: Dense = rng.unitary2().0.iter().map(|r| r.to_vec()).collect();
            let u4: Dense = rng.unitary4().0.iter().map(|r| r.to_vec()).collect();
            for u in [u2, u4] {
                let d = decompose_unitary(&u, &native(&["rz", "sx", "x", "cx"])).unwrap();
                assert!(d.error < 1e-9);
                assert!(
                    d.instructions
                        .iter()
                        .filter(|i| i.qubits.len() == 2)
                        .count()
                        <= 3
                );
            }
        }
    }

    #[test]
    fn errors_are_explained() {
        let set = native(&["rz", "sx", "cx"]);
        assert!(
            decompose_gate("nope", &[], &set)
                .unwrap_err()
                .contains("unknown gate")
        );
        assert!(
            decompose_gate("rx", &[], &set)
                .unwrap_err()
                .contains("the gate rx takes 1 parameter")
        );
        let not_unitary = vec![vec![ONE, ONE], vec![ZERO, ONE]];
        assert!(
            decompose_unitary(&not_unitary, &set)
                .unwrap_err()
                .contains("not unitary")
        );
        assert!(decompose_gate("cx", &[], &native(&["rz", "sx", "swap"])).is_err());
        let identity: Dense = (0..4)
            .map(|i| (0..4).map(|j| if i == j { ONE } else { ZERO }).collect())
            .collect();
        assert!(
            decompose_unitary(&identity, &native(&["rz", "sx"]))
                .unwrap_err()
                .contains("no supported two-qubit gate")
        );
    }
}
