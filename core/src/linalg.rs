//! Fixed-size complex matrices: 2×2 for single-qubit gates, 4×4 for
//! two-qubit gates, with a Jacobi eigensolver for real symmetric 4×4 matrices.
//!
//! Two-qubit matrices are big-endian: in `kron(a, b)`, `a` acts on the first
//! (most significant) qubit.

use num_complex::Complex64;
use std::ops::{Mul, MulAssign};

pub type C = Complex64;

pub const ZERO: C = C::new(0.0, 0.0);
pub const ONE: C = C::new(1.0, 0.0);
pub const I: C = C::new(0.0, 1.0);

/// A 2×2 complex matrix.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct M2(pub [[C; 2]; 2]);

/// A 4×4 complex matrix.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct M4(pub [[C; 4]; 4]);

impl M2 {
    pub const IDENTITY: M2 = M2([[ONE, ZERO], [ZERO, ONE]]);

    pub fn new(a: C, b: C, c: C, d: C) -> M2 {
        M2([[a, b], [c, d]])
    }

    pub fn adjoint(&self) -> M2 {
        let m = &self.0;
        M2([
            [m[0][0].conj(), m[1][0].conj()],
            [m[0][1].conj(), m[1][1].conj()],
        ])
    }

    pub fn det(&self) -> C {
        let m = &self.0;
        m[0][0] * m[1][1] - m[0][1] * m[1][0]
    }

    pub fn scale(&self, s: C) -> M2 {
        let m = &self.0;
        M2([[m[0][0] * s, m[0][1] * s], [m[1][0] * s, m[1][1] * s]])
    }

    /// Whether `self` and `other` are equal element by element within `tol`.
    pub fn close(&self, other: &M2, tol: f64) -> bool {
        (0..2).all(|i| (0..2).all(|j| (self.0[i][j] - other.0[i][j]).norm() <= tol))
    }

    /// Whether `self` is the identity up to a global phase.
    pub fn is_identity_up_to_phase(&self, tol: f64) -> bool {
        let m = &self.0;
        m[0][1].norm() < tol && m[1][0].norm() < tol && (m[0][0] - m[1][1]).norm() < tol
    }
}

impl Mul for M2 {
    type Output = M2;
    fn mul(self, rhs: M2) -> M2 {
        let (a, b) = (&self.0, &rhs.0);
        let mut out = [[ZERO; 2]; 2];
        for (i, row) in out.iter_mut().enumerate() {
            for (j, value) in row.iter_mut().enumerate() {
                *value = a[i][0] * b[0][j] + a[i][1] * b[1][j];
            }
        }
        M2(out)
    }
}

impl M4 {
    pub const IDENTITY: M4 = M4([
        [ONE, ZERO, ZERO, ZERO],
        [ZERO, ONE, ZERO, ZERO],
        [ZERO, ZERO, ONE, ZERO],
        [ZERO, ZERO, ZERO, ONE],
    ]);

    pub fn from_real(m: &[[f64; 4]; 4]) -> M4 {
        let mut out = [[ZERO; 4]; 4];
        for i in 0..4 {
            for j in 0..4 {
                out[i][j] = C::new(m[i][j], 0.0);
            }
        }
        M4(out)
    }

    pub fn adjoint(&self) -> M4 {
        let mut out = [[ZERO; 4]; 4];
        for i in 0..4 {
            for j in 0..4 {
                out[i][j] = self.0[j][i].conj();
            }
        }
        M4(out)
    }

    pub fn transpose(&self) -> M4 {
        let mut out = [[ZERO; 4]; 4];
        for i in 0..4 {
            for j in 0..4 {
                out[i][j] = self.0[j][i];
            }
        }
        M4(out)
    }

    pub fn scale(&self, s: C) -> M4 {
        let mut out = self.0;
        for row in out.iter_mut() {
            for value in row.iter_mut() {
                *value *= s;
            }
        }
        M4(out)
    }

    /// The determinant, by Gaussian elimination with partial pivoting.
    pub fn det(&self) -> C {
        let mut m = self.0;
        let mut det = ONE;
        for col in 0..4 {
            let pivot = (col..4)
                .max_by(|&a, &b| m[a][col].norm().total_cmp(&m[b][col].norm()))
                .unwrap();
            if m[pivot][col].norm() == 0.0 {
                return ZERO;
            }
            if pivot != col {
                m.swap(pivot, col);
                det = -det;
            }
            det *= m[col][col];
            for row in col + 1..4 {
                let factor = m[row][col] / m[col][col];
                for k in col..4 {
                    let value = m[col][k];
                    m[row][k] -= factor * value;
                }
            }
        }
        det
    }

    /// `self` with its two qubits exchanged: SWAP · self · SWAP.
    pub fn swapped(&self) -> M4 {
        const P: [usize; 4] = [0, 2, 1, 3];
        let mut out = [[ZERO; 4]; 4];
        for i in 0..4 {
            for j in 0..4 {
                out[i][j] = self.0[P[i]][P[j]];
            }
        }
        M4(out)
    }

    pub fn close(&self, other: &M4, tol: f64) -> bool {
        (0..4).all(|i| (0..4).all(|j| (self.0[i][j] - other.0[i][j]).norm() <= tol))
    }

    pub fn trace(&self) -> C {
        (0..4).map(|i| self.0[i][i]).sum()
    }

    /// Frobenius inner product: trace(self† · other).
    pub fn inner(&self, other: &M4) -> C {
        let mut sum = ZERO;
        for i in 0..4 {
            for j in 0..4 {
                sum += self.0[i][j].conj() * other.0[i][j];
            }
        }
        sum
    }
}

impl Mul for M4 {
    type Output = M4;
    fn mul(self, rhs: M4) -> M4 {
        let (a, b) = (&self.0, &rhs.0);
        let mut out = [[ZERO; 4]; 4];
        for i in 0..4 {
            for j in 0..4 {
                out[i][j] =
                    a[i][0] * b[0][j] + a[i][1] * b[1][j] + a[i][2] * b[2][j] + a[i][3] * b[3][j];
            }
        }
        M4(out)
    }
}

impl MulAssign for M4 {
    fn mul_assign(&mut self, rhs: M4) {
        *self = *self * rhs;
    }
}

/// `a ⊗ b`, with `a` on the first (most significant) qubit.
pub fn kron(a: &M2, b: &M2) -> M4 {
    let mut out = [[ZERO; 4]; 4];
    for i in 0..2 {
        for j in 0..2 {
            for k in 0..2 {
                for l in 0..2 {
                    out[2 * i + k][2 * j + l] = a.0[i][j] * b.0[k][l];
                }
            }
        }
    }
    M4(out)
}

/// Whether `a = e^{iφ} b` for some φ, within `tol` per element.
pub fn equal_up_to_phase4(a: &M4, b: &M4, tol: f64) -> bool {
    let overlap = b.inner(a);
    if overlap.norm() < 1e-12 {
        return false;
    }
    let phase = overlap / overlap.norm();
    a.close(&b.scale(phase), tol)
}

/// Eigenvectors of a real symmetric 4×4 matrix by the cyclic Jacobi method,
/// as the columns of an orthogonal matrix.
pub fn symmetric_eigenvectors(m: &[[f64; 4]; 4]) -> [[f64; 4]; 4] {
    let mut a = *m;
    let mut v = [[0.0; 4]; 4];
    for (i, row) in v.iter_mut().enumerate() {
        row[i] = 1.0;
    }
    for _ in 0..100 {
        let off: f64 = (0..4)
            .flat_map(|i| (0..4).filter(move |&j| j != i).map(move |j| (i, j)))
            .map(|(i, j)| a[i][j] * a[i][j])
            .sum();
        if off < 1e-30 {
            break;
        }
        for p in 0..3 {
            for q in p + 1..4 {
                if a[p][q].abs() < 1e-300 {
                    continue;
                }
                let theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let t = if theta == 0.0 { 1.0 } else { t };
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                for k in 0..4 {
                    let (akp, akq) = (a[k][p], a[k][q]);
                    a[k][p] = c * akp - s * akq;
                    a[k][q] = s * akp + c * akq;
                }
                for k in 0..4 {
                    let (apk, aqk) = (a[p][k], a[q][k]);
                    a[p][k] = c * apk - s * aqk;
                    a[q][k] = s * apk + c * aqk;
                }
                for row in v.iter_mut() {
                    let (vkp, vkq) = (row[p], row[q]);
                    row[p] = c * vkp - s * vkq;
                    row[q] = s * vkp + c * vkq;
                }
            }
        }
    }
    v
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn jacobi_diagonalises() {
        let m = [
            [4.0, 1.0, -2.0, 0.5],
            [1.0, 3.0, 0.0, 1.5],
            [-2.0, 0.0, 1.0, -1.0],
            [0.5, 1.5, -1.0, 2.0],
        ];
        let v = symmetric_eigenvectors(&m);
        // vᵀ m v is diagonal and vᵀ v = 1.
        for i in 0..4 {
            for j in 0..4 {
                let d: f64 = (0..4)
                    .flat_map(|k| (0..4).map(move |l| (k, l)))
                    .map(|(k, l)| v[k][i] * m[k][l] * v[l][j])
                    .sum();
                let o: f64 = (0..4).map(|k| v[k][i] * v[k][j]).sum();
                if i != j {
                    assert!(d.abs() < 1e-12, "{d}");
                    assert!(o.abs() < 1e-12);
                } else {
                    assert!((o - 1.0).abs() < 1e-12);
                }
            }
        }
    }

    #[test]
    fn determinant() {
        let a = M2::new(C::new(0.0, 1.0), ONE, ZERO, C::new(2.0, 0.0));
        let b = M2::new(ONE, C::new(1.0, 1.0), C::new(0.5, 0.0), C::new(-1.0, 0.0));
        let k = kron(&a, &b);
        // det(a ⊗ b) = det(a)² det(b)².
        let expected = a.det() * a.det() * b.det() * b.det();
        assert!((k.det() - expected).norm() < 1e-12);
    }
}
