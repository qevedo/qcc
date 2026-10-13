//! Decompositions into any entangling two-qubit gate, found numerically.
//!
//! For gates without a closed-form template (iSWAP, √iSWAP, √SWAP, …) this
//! looks, for k = 1, 2, 3, … uses of the gate `G`, for single-qubit layers
//! with `(A₁⊗A₀) G L_{k−1} G ⋯ L₁ G (B₁⊗B₀) = e^{iφ} U`, in two stages:
//!
//! 1. Two gates are equal up to single-qubit gates exactly when their Makhlin
//!    invariants (G1 complex, G2 real) agree. Solving those three equations
//!    for the inner layers is cheap and decides whether k uses can work.
//! 2. Near targets with repeated Weyl coordinates (controlled rotations, for
//!    example) the invariants pin the point down only to about the square
//!    root of their error, so the decomposition is polished against the
//!    target matrix itself: all layers, the outer gates (started from the KAK
//!    decompositions) and the phase, until it is exact to rounding.
//!
//! Both stages use Levenberg–Marquardt from seeded random starts. The first k
//! that works is used, so the count is minimal unless every start for a
//! smaller k missed.

use crate::gates::{ry, rz};
use crate::linalg::{C, M2, M4, kron};
use crate::one_qubit::euler_zyz;
use crate::rng::Rng;
use crate::two_qubit::{Kak, Step, canonical, kak, magic};

/// The most uses of the gate the search tries.
pub const MAX_USES: usize = 6;
const STARTS: usize = 40;

/// Makhlin invariants (Re G1, Im G1, G2) of a two-qubit unitary.
pub fn invariants(u: &M4) -> [f64; 3] {
    let (q, q_dag) = magic();
    let ub = q_dag * *u * q;
    let m = ub.transpose() * ub;
    let det = u.det();
    let tr = m.trace();
    let tr2 = (m * m).trace();
    let g1 = tr * tr / (det * 16.0);
    let g2 = (tr * tr - tr2) / (det * 4.0);
    [g1.re, g1.im, g2.re]
}

fn euler(p: &[f64]) -> M2 {
    rz(p[0]) * ry(p[1]) * rz(p[2])
}

fn euler_params(m: &M2) -> [f64; 3] {
    let (theta, phi, lam) = euler_zyz(m);
    [phi, theta, lam]
}

/// G L_{k−1} G ⋯ L₁ G for the inner layers in `params` (6 numbers each).
fn inner(g: &M4, params: &[f64]) -> M4 {
    let mut out = *g;
    for chunk in params.chunks(6) {
        out = *g * kron(&euler(&chunk[0..3]), &euler(&chunk[3..6])) * out;
    }
    out
}

/// Levenberg–Marquardt on `residual(x)`, from `x`. Returns the solution and
/// the squared norm of its residual.
fn least_squares(
    residual: impl Fn(&[f64]) -> Vec<f64>,
    mut x: Vec<f64>,
    iterations: usize,
) -> (Vec<f64>, f64) {
    let n = x.len();
    let norm2 = |r: &[f64]| r.iter().map(|v| v * v).sum::<f64>();
    let mut r = residual(&x);
    let mut f = norm2(&r);
    let mut mu = 1e-3;
    for _ in 0..iterations {
        if f < 1e-30 {
            break;
        }
        let m = r.len();
        // Central-difference Jacobian, column by column.
        let mut jac = vec![vec![0.0; m]; n];
        for (k, column) in jac.iter_mut().enumerate() {
            let h = 1e-7;
            let (mut plus, mut minus) = (x.clone(), x.clone());
            plus[k] += h;
            minus[k] -= h;
            let (rp, rm) = (residual(&plus), residual(&minus));
            for i in 0..m {
                column[i] = (rp[i] - rm[i]) / (2.0 * h);
            }
        }
        let jtj: Vec<Vec<f64>> = (0..n)
            .map(|a| {
                (0..n)
                    .map(|b| (0..m).map(|i| jac[a][i] * jac[b][i]).sum())
                    .collect()
            })
            .collect();
        let jtr: Vec<f64> = (0..n)
            .map(|a| (0..m).map(|i| jac[a][i] * r[i]).sum())
            .collect();
        let mut improved = false;
        for _ in 0..16 {
            let mut a = jtj.clone();
            for (i, row) in a.iter_mut().enumerate() {
                row[i] += mu * (1.0 + jtj[i][i]);
            }
            let Some(step) = solve(a, jtr.iter().map(|v| -v).collect()) else {
                break;
            };
            let candidate: Vec<f64> = x.iter().zip(&step).map(|(a, b)| a + b).collect();
            let rc = residual(&candidate);
            let fc = norm2(&rc);
            if fc < f {
                (x, r, f) = (candidate, rc, fc);
                mu = (mu / 4.0).max(1e-18);
                improved = true;
                break;
            }
            mu *= 8.0;
        }
        if !improved {
            break;
        }
    }
    (x, f)
}

/// Solve `a x = b` by Gaussian elimination with partial pivoting.
fn solve(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Option<Vec<f64>> {
    let n = b.len();
    for col in 0..n {
        let pivot = (col..n).max_by(|&i, &j| a[i][col].abs().total_cmp(&a[j][col].abs()))?;
        if a[pivot][col].abs() < 1e-300 {
            return None;
        }
        a.swap(pivot, col);
        b.swap(pivot, col);
        for row in col + 1..n {
            let f = a[row][col] / a[col][col];
            if f == 0.0 {
                continue;
            }
            for k in col..n {
                a[row][k] -= f * a[col][k];
            }
            b[row] -= f * b[col];
        }
    }
    let mut x = vec![0.0; n];
    for row in (0..n).rev() {
        let s: f64 = (row + 1..n).map(|k| a[row][k] * x[k]).sum();
        x[row] = (b[row] - s) / a[row][row];
    }
    Some(x)
}

/// The fewest uses of gate `g` (named `name`) that build `u`, as steps (the
/// outer single-qubit gates included); `None` when no decomposition with up
/// to `MAX_USES` uses was found.
pub fn numeric_decomposition(u: &M4, target: &Kak, name: &str, g: &M4) -> Option<Vec<Step>> {
    let point = target.point;
    if point.iter().all(|c| c.abs() < 1e-9) {
        // A product of single-qubit gates.
        return Some(vec![
            Step::U(0, target.a1 * target.b1),
            Step::U(1, target.a0 * target.b0),
        ]);
    }
    let goal = invariants(&canonical(point));
    // Seeded from the target, so a decomposition is reproducible.
    let seed = point.iter().fold(0x5EED_u64, |s, c| {
        s.rotate_left(17) ^ ((c * 1e9).round() as i64 as u64)
    });
    let mut rng = Rng::new(seed);
    for uses in 1..=MAX_USES {
        let n = 6 * (uses - 1);
        let starts = if uses == 1 { 1 } else { STARTS };
        for attempt in 0..starts {
            // Identity layers first: they find exact powers (√SWAP² = SWAP) with
            // no single-qubit gates in between.
            let start: Vec<f64> = if attempt == 0 {
                vec![0.0; n]
            } else {
                (0..n)
                    .map(|_| rng.range(-std::f64::consts::PI, std::f64::consts::PI))
                    .collect()
            };
            let residual = |x: &[f64]| {
                let inv = invariants(&inner(g, x));
                vec![inv[0] - goal[0], inv[1] - goal[1], inv[2] - goal[2]]
            };
            let (x, f) = least_squares(residual, start, 300);
            if f > 1e-8 {
                continue;
            }
            if let Some(steps) = polish(u, g, name, &x) {
                return Some(steps);
            }
        }
    }
    None
}

/// Refine inner layers `x` that nearly reach `u`'s class into an exact
/// decomposition of `u`, outer gates and phase included.
fn polish(u: &M4, g: &M4, name: &str, x: &[f64]) -> Option<Vec<Step>> {
    // Outer gates from the two KAK decompositions, as in `two_qubit::instantiate`.
    let (target, reference) = (kak(u), kak(&inner(g, x)));
    let before = [
        reference.b1.adjoint() * target.b1,
        reference.b0.adjoint() * target.b0,
    ];
    let after = [
        target.a1 * reference.a1.adjoint(),
        target.a0 * reference.a0.adjoint(),
    ];
    let n = x.len();
    let mut start: Vec<f64> = x.to_vec();
    for m in before.iter().chain(&after) {
        start.extend(euler_params(m));
    }
    start.push(0.0);
    let full = |p: &[f64]| -> M4 {
        let b = kron(&euler(&p[n..n + 3]), &euler(&p[n + 3..n + 6]));
        let a = kron(&euler(&p[n + 6..n + 9]), &euler(&p[n + 9..n + 12]));
        (a * inner(g, &p[..n]) * b).scale(C::from_polar(1.0, p[n + 12]))
    };
    // Start the phase from the overlap with u.
    start[n + 12] = -u.inner(&full(&start)).arg();
    let residual = |p: &[f64]| {
        let f = full(p);
        let mut r = Vec::with_capacity(32);
        for i in 0..4 {
            for j in 0..4 {
                let d = f.0[i][j] - u.0[i][j];
                r.push(d.re);
                r.push(d.im);
            }
        }
        r
    };
    let (p, f) = least_squares(residual, start, 100);
    if f > 1e-26 {
        return None;
    }
    let gate = || Step::Gate(name.to_string(), vec![], [0, 1]);
    let mut steps = vec![
        Step::U(0, euler(&p[n..n + 3])),
        Step::U(1, euler(&p[n + 3..n + 6])),
        gate(),
    ];
    for chunk in p[..n].chunks(6) {
        steps.extend([
            Step::U(0, euler(&chunk[0..3])),
            Step::U(1, euler(&chunk[3..6])),
            gate(),
        ]);
    }
    steps.extend([
        Step::U(0, euler(&p[n + 6..n + 9])),
        Step::U(1, euler(&p[n + 9..n + 12])),
    ]);
    Some(steps)
}

/// Whether `g` can build every two-qubit gate with single-qubit gates: any
/// entangling gate except those equal to SWAP up to single-qubit gates.
pub fn is_universal(point: [f64; 3]) -> bool {
    let q = std::f64::consts::FRAC_PI_4;
    let zero = point.iter().all(|c| c.abs() < 1e-9);
    let swap = point.iter().all(|c| (c - q).abs() < 1e-9);
    !zero && !swap
}
