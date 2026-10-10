//! Helpers for the unit tests: a seeded generator and random unitaries.

use crate::linalg::{C, M2, M4, ZERO};

/// SplitMix64, enough for test inputs.
pub struct Rng(u64);

impl Rng {
    pub fn new(seed: u64) -> Rng {
        Rng(seed)
    }

    pub fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Uniform in [0, 1).
    pub fn uniform(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    pub fn range(&mut self, low: f64, high: f64) -> f64 {
        low + (high - low) * self.uniform()
    }

    pub fn below(&mut self, n: usize) -> usize {
        (self.next_u64() % n as u64) as usize
    }

    fn normal(&mut self) -> f64 {
        let (u1, u2) = (self.uniform().max(1e-300), self.uniform());
        (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos()
    }

    fn columns<const N: usize>(&mut self) -> [[C; N]; N] {
        // Gram–Schmidt on Gaussian columns gives a Haar-random unitary.
        let mut cols = [[ZERO; N]; N];
        for j in 0..N {
            let mut v = [ZERO; N];
            for value in v.iter_mut() {
                *value = C::new(self.normal(), self.normal());
            }
            for prev in cols.iter().take(j) {
                let dot: C = (0..N).map(|i| prev[i].conj() * v[i]).sum();
                for i in 0..N {
                    v[i] -= dot * prev[i];
                }
            }
            let norm = v.iter().map(|x| x.norm_sqr()).sum::<f64>().sqrt();
            for i in 0..N {
                cols[j][i] = v[i] / norm;
            }
        }
        cols
    }

    pub fn unitary2(&mut self) -> M2 {
        let c = self.columns::<2>();
        M2([[c[0][0], c[1][0]], [c[0][1], c[1][1]]])
    }

    pub fn unitary4(&mut self) -> M4 {
        let c = self.columns::<4>();
        let mut m = [[ZERO; 4]; 4];
        for i in 0..4 {
            for j in 0..4 {
                m[i][j] = c[j][i];
            }
        }
        M4(m)
    }
}

pub fn equal_up_to_phase2(a: &M2, b: &M2, tol: f64) -> bool {
    let overlap: C = (0..2)
        .flat_map(|i| (0..2).map(move |j| (i, j)))
        .map(|(i, j)| b.0[i][j].conj() * a.0[i][j])
        .sum();
    if overlap.norm() < 1e-12 {
        return false;
    }
    a.close(&b.scale(overlap / overlap.norm()), tol)
}

/// The dense unitary of `instructions` on `n` qubits, qubit 0 most significant.
pub fn circuit_unitary(instructions: &[crate::ir::Instruction], n: usize) -> Vec<Vec<C>> {
    use crate::gates::{matrix_1q, matrix_2q};
    let dim = 1 << n;
    let mut u: Vec<Vec<C>> = (0..dim)
        .map(|i| {
            (0..dim)
                .map(|j| if i == j { C::new(1.0, 0.0) } else { ZERO })
                .collect()
        })
        .collect();
    let bit = |index: usize, q: usize| (index >> (n - 1 - q)) & 1;
    for inst in instructions {
        if inst.name == "barrier" {
            continue;
        }
        let k = inst.qubits.len();
        // Gate matrix on its own qubits, big-endian.
        let g: Vec<Vec<C>> = if k == 1 {
            matrix_1q(&inst.name, &inst.params)
                .unwrap()
                .0
                .iter()
                .map(|r| r.to_vec())
                .collect()
        } else {
            matrix_2q(&inst.name, &inst.params)
                .unwrap()
                .0
                .iter()
                .map(|r| r.to_vec())
                .collect()
        };
        let mut next = vec![vec![ZERO; dim]; dim];
        for row in 0..dim {
            let sub_row: usize = inst
                .qubits
                .iter()
                .fold(0, |acc, &q| (acc << 1) | bit(row, q));
            for sub_col in 0..(1 << k) {
                let coefficient = g[sub_row][sub_col];
                if coefficient == ZERO {
                    continue;
                }
                let mut source = row;
                for (pos, &q) in inst.qubits.iter().enumerate() {
                    let value = (sub_col >> (k - 1 - pos)) & 1;
                    let mask = 1 << (n - 1 - q);
                    source = if value == 1 {
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

pub fn dense_equal_up_to_phase(a: &[Vec<C>], b: &[Vec<C>], tol: f64) -> bool {
    let overlap: C = a
        .iter()
        .zip(b)
        .flat_map(|(ra, rb)| ra.iter().zip(rb).map(|(x, y)| y.conj() * x))
        .sum();
    if overlap.norm() < 1e-12 {
        return false;
    }
    let phase = overlap / overlap.norm();
    a.iter().zip(b).all(|(ra, rb)| {
        ra.iter()
            .zip(rb)
            .all(|(x, y)| (x - y * phase).norm() <= tol)
    })
}
