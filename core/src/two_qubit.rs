//! Optimal two-qubit synthesis (KAK / Cartan decomposition).
//!
//! Every two-qubit unitary is `e^{iγ} (A1 ⊗ A0) · exp(i(a XX + b YY + c ZZ)) ·
//! (B1 ⊗ B0)` with a point `(a, b, c)` of the Weyl chamber `π/4 ≥ a ≥ b ≥ |c|`.
//! The point fixes the minimal number of native two-qubit gates: 0–3 for
//! CX-like gates, one per non-zero coordinate for Ising rotations. A template
//! circuit with the same point is dressed with single-qubit gates from the two
//! KAK decompositions. Matches `qevedo/compiler/synthesis/two_qubit.py`.

use crate::gates::{cx, matrix_2q, pauli, rotation, rx, ry, rz, sx};
use crate::linalg::{C, I, M2, M4, ONE, ZERO, equal_up_to_phase4, kron, symmetric_eigenvectors};
use crate::one_qubit::{euler_zyz, remainder};
use std::collections::HashMap;
use std::f64::consts::{FRAC_1_SQRT_2, FRAC_PI_2, FRAC_PI_4, PI, TAU};
use std::sync::LazyLock;

pub const TOL: f64 = 1e-9;

/// `numpy.allclose` for one value: |a − b| ≤ atol + 1e-5 |b|.
fn close(a: f64, b: f64, atol: f64) -> bool {
    (a - b).abs() <= atol + 1e-5 * b.abs()
}

fn close3(a: [f64; 3], b: [f64; 3], atol: f64) -> bool {
    (0..3).all(|i| close(a[i], b[i], atol))
}

static MAGIC: LazyLock<M4> = LazyLock::new(|| {
    let r = C::new(FRAC_1_SQRT_2, 0.0);
    let (o, i) = (ONE * r, I * r);
    M4([
        [o, ZERO, ZERO, i],
        [ZERO, i, o, ZERO],
        [ZERO, i, -o, ZERO],
        [o, ZERO, ZERO, -i],
    ])
});
static MAGIC_DAG: LazyLock<M4> = LazyLock::new(|| MAGIC.adjoint());
/// The diagonals of XX, YY and ZZ in the magic basis (each ±1).
static DIAG: LazyLock<[[f64; 4]; 3]> = LazyLock::new(|| {
    let mut out = [[0.0; 4]; 3];
    for (axis, row) in out.iter_mut().enumerate() {
        let p = pauli(axis);
        let d = *MAGIC_DAG * kron(&p, &p) * *MAGIC;
        for (k, value) in row.iter_mut().enumerate() {
            *value = d.0[k][k].re;
        }
    }
    out
});

/// exp(i(a XX + b YY + c ZZ)).
pub fn canonical(point: [f64; 3]) -> M4 {
    let mut d = M4([[ZERO; 4]; 4]);
    for k in 0..4 {
        let angle: f64 = (0..3).map(|axis| DIAG[axis][k] * point[axis]).sum();
        d.0[k][k] = C::from_polar(1.0, angle);
    }
    *MAGIC * d * *MAGIC_DAG
}

/// `U = e^{i phase} (a1 ⊗ a0) canonical(point) (b1 ⊗ b0)`; `a1` and `b1` act on
/// the first (most significant) qubit.
#[derive(Clone, Debug)]
pub struct Kak {
    pub phase: f64,
    pub point: [f64; 3],
    pub a1: M2,
    pub a0: M2,
    pub b1: M2,
    pub b0: M2,
}

impl Kak {
    pub fn matrix(&self) -> M4 {
        (kron(&self.a1, &self.a0) * canonical(self.point) * kron(&self.b1, &self.b0))
            .scale(C::from_polar(1.0, self.phase))
    }
}

fn su2(m: M2) -> M2 {
    m.scale(m.det().sqrt().inv())
}

/// Split `k = phase · (p ⊗ q)` with `p` and `q` in SU(2).
fn factor(k: &M4) -> (M2, M2, C) {
    // m[(a, a'), (c, c')] = k[(a, c), (a', c')] is the rank-one vec(p) vec(q)ᵀ.
    let entry = |r: usize, c: usize| k.0[2 * (r / 2) + c / 2][2 * (r % 2) + c % 2];
    let (mut r0, mut c0) = (0, 0);
    for r in 0..4 {
        for c in 0..4 {
            if entry(r, c).norm() > entry(r0, c0).norm() {
                (r0, c0) = (r, c);
            }
        }
    }
    let p = su2(M2([
        [entry(0, c0), entry(1, c0)],
        [entry(2, c0), entry(3, c0)],
    ]));
    let q = su2(M2([
        [entry(r0, 0), entry(r0, 1)],
        [entry(r0, 2), entry(r0, 3)],
    ]));
    let scale = kron(&p, &q).inner(k) / 4.0;
    (p, q, scale)
}

/// A real orthogonal `P` with det 1 that diagonalises the symmetric unitary `m`.
///
/// The real and imaginary parts of `m` commute, so a generic real combination
/// of them has the common eigenvectors.
fn orthogonal_eigenvectors(m: &M4) -> [[f64; 4]; 4] {
    let mut r = 0.5f64;
    for _ in 0..100 {
        let mut combo = [[0.0; 4]; 4];
        for i in 0..4 {
            for j in 0..4 {
                combo[i][j] = r * m.0[i][j].re + (1.0 - r) * m.0[i][j].im;
            }
        }
        let mut p = symmetric_eigenvectors(&combo);
        let pm = M4::from_real(&p);
        let d = pm.transpose() * *m * pm;
        let diagonal = (0..4).all(|i| (0..4).all(|j| i == j || d.0[i][j].norm() <= 1e-10));
        if diagonal {
            if pm.det().re < 0.0 {
                for row in p.iter_mut() {
                    row[0] = -row[0];
                }
            }
            return p;
        }
        // The golden-ratio sequence covers (0, 1) evenly.
        r = (r + 0.618_033_988_749_895) % 1.0;
    }
    panic!("could not diagonalise the KAK matrix");
}

/// The KAK decomposition of a 4×4 unitary, with its point in the Weyl chamber.
pub fn kak(u: &M4) -> Kak {
    let g = u.det().arg() / 4.0;
    let up = *MAGIC_DAG * u.scale(C::from_polar(1.0, -g)) * *MAGIC;
    let m = up.transpose() * up;
    let p = orthogonal_eigenvectors(&m);
    let pm = M4::from_real(&p);
    let dm = pm.transpose() * m * pm;
    let mut d: [C; 4] = std::array::from_fn(|k| dm.0[k][k].sqrt());
    let mut conj_d = M4([[ZERO; 4]; 4]);
    for k in 0..4 {
        conj_d.0[k][k] = d[k].conj();
    }
    let mut k1 = up * pm * conj_d;
    if k1.det().re < 0.0 {
        d[0] = -d[0];
        for row in k1.0.iter_mut() {
            row[0] = -row[0];
        }
    }
    // up = k1 · diag(d) · pᵀ with k1 and pᵀ in SO(4).
    let theta: [f64; 4] = std::array::from_fn(|k| d[k].arg());
    let offset = theta.iter().sum::<f64>() / 4.0;
    let point: [f64; 3] = std::array::from_fn(|axis| {
        (0..4)
            .map(|k| DIAG[axis][k] * (theta[k] - offset))
            .sum::<f64>()
            / 4.0
    });
    let mut k1_real = k1;
    for row in k1_real.0.iter_mut() {
        for value in row.iter_mut() {
            *value = C::new(value.re, 0.0);
        }
    }
    let (left1, left0, s1) = factor(&(*MAGIC * k1_real * *MAGIC_DAG));
    let (right1, right0, s2) = factor(&(*MAGIC * pm.transpose() * *MAGIC_DAG));
    let phase = g + offset + s1.arg() + s2.arg();
    canonicalize(phase, point, [left1, left0], [right1, right0])
}

/// exp(-iπ/4 P) for the Pauli on `axis`.
fn quarter_turn(axis: usize) -> M2 {
    let p = pauli(axis);
    let r = FRAC_1_SQRT_2;
    let mut out = p.scale(C::new(0.0, -r));
    out.0[0][0] += r;
    out.0[1][1] += r;
    out
}

/// Move `v` into the Weyl chamber, keeping `U` unchanged.
fn canonicalize(mut phase: f64, mut v: [f64; 3], mut left: [M2; 2], mut right: [M2; 2]) -> Kak {
    // exp(i(v + step·π/2) PP) = exp(i v PP) · (i PP)^step
    let shift = |v: &mut [f64; 3], right: &mut [M2; 2], phase: &mut f64, k: usize, step: f64| {
        v[k] += step * FRAC_PI_2;
        right[0] = pauli(k) * right[0];
        right[1] = pauli(k) * right[1];
        *phase -= step * FRAC_PI_2;
    };
    // (P ⊗ I) canonical(v) (P ⊗ I) flips the two coordinates P anticommutes with.
    let negate =
        |v: &mut [f64; 3], left: &mut [M2; 2], right: &mut [M2; 2], k1: usize, k2: usize| {
            let p = pauli(3 - k1 - k2);
            v[k1] = -v[k1];
            v[k2] = -v[k2];
            left[0] = left[0] * p;
            right[0] = p * right[0];
        };
    // A quarter turn about the third axis on both qubits exchanges the other two.
    let swap = |v: &mut [f64; 3], left: &mut [M2; 2], right: &mut [M2; 2], k1: usize, k2: usize| {
        let r = quarter_turn(3 - k1 - k2);
        v.swap(k1, k2);
        left[0] = left[0] * r.adjoint();
        left[1] = left[1] * r.adjoint();
        right[0] = r * right[0];
        right[1] = r * right[1];
    };
    for k in 0..3 {
        while v[k] > FRAC_PI_4 + TOL {
            shift(&mut v, &mut right, &mut phase, k, -1.0);
        }
        while v[k] <= -FRAC_PI_4 + TOL {
            shift(&mut v, &mut right, &mut phase, k, 1.0);
        }
    }
    for _ in 0..3 {
        for k in 0..2 {
            if v[k].abs() < v[k + 1].abs() - TOL {
                swap(&mut v, &mut left, &mut right, k, k + 1);
            }
        }
    }
    if v[0] < -TOL {
        negate(&mut v, &mut left, &mut right, 0, 2);
    }
    if v[1] < -TOL {
        negate(&mut v, &mut left, &mut right, 1, 2);
    }
    if (v[0] - FRAC_PI_4).abs() < TOL && v[2] < -TOL {
        shift(&mut v, &mut right, &mut phase, 0, -1.0);
        negate(&mut v, &mut left, &mut right, 0, 2);
    }
    Kak {
        phase,
        point: v,
        a1: left[0],
        a0: left[1],
        b1: right[0],
        b0: right[1],
    }
}

// -- synthesis ---------------------------------------------------------------------

/// A step of a two-qubit sequence, in time order. Qubit 0 is the first operand
/// of the decomposed gate.
#[derive(Clone, Debug)]
pub enum Step {
    /// A single-qubit unitary still to be synthesized.
    U(usize, M2),
    /// A native two-qubit gate on `[first, second]`.
    Gate(String, Vec<f64>, [usize; 2]),
}

/// A native two-qubit gate and how the templates use it.
#[derive(Clone, Debug)]
pub struct TwoQubitBasis {
    pub name: String,
    pub kind: BasisKind,
}

#[derive(Clone, Debug)]
pub enum BasisKind {
    /// Equal to CX up to single-qubit gates: `CX ∝ (after1 ⊗ after0) · G · (before1 ⊗ before0)`,
    /// stored as `[after1, after0, before1, before0]`.
    CxLike(Box<[M2; 4]>),
    /// exp(-iθ/2 P⊗Q) with the Paulis `(first, second)`.
    Ising(usize, usize),
}

/// How to use native gate `name`, or None when no template uses it yet.
pub fn two_qubit_basis(name: &str) -> Option<TwoQubitBasis> {
    let ising = match name {
        "rxx" => Some((0, 0)),
        "ryy" => Some((1, 1)),
        "rzz" => Some((2, 2)),
        "rzx" => Some((2, 0)),
        _ => None,
    };
    if let Some((p, q)) = ising {
        return Some(TwoQubitBasis {
            name: name.to_string(),
            kind: BasisKind::Ising(p, q),
        });
    }
    let g = matrix_2q(name, &[])?;
    let kg = kak(&g);
    if !close3(kg.point, [FRAC_PI_4, 0.0, 0.0], TOL) {
        return None;
    }
    let kc = kak(&cx());
    // CX ∝ (Ac Ag†) G (Bg† Bc) on each qubit.
    let locals = [
        kc.a1 * kg.a1.adjoint(),
        kc.a0 * kg.a0.adjoint(),
        kg.b1.adjoint() * kc.b1,
        kg.b0.adjoint() * kc.b0,
    ];
    Some(TwoQubitBasis {
        name: name.to_string(),
        kind: BasisKind::CxLike(Box::new(locals)),
    })
}

/// The native two-qubit gates some template can use, by name.
pub fn bases_2q<'a>(native: impl IntoIterator<Item = &'a String>) -> Vec<TwoQubitBasis> {
    let mut names: Vec<&String> = native.into_iter().collect();
    names.sort();
    names
        .into_iter()
        .filter_map(|n| two_qubit_basis(n))
        .collect()
}

fn cx_step(control: usize, target: usize) -> Step {
    Step::Gate("cx".to_string(), vec![], [control, target])
}

/// Circuits with the fewest CXs whose Weyl point is `point`, in variants that
/// differ by rotations commuting with a neighbouring CX.
fn cx_templates(point: [f64; 3]) -> Vec<Vec<Step>> {
    let [a, b, c] = point;
    if close3(point, [0.0; 3], TOL) {
        return vec![vec![]];
    }
    if close3(point, [FRAC_PI_4, 0.0, 0.0], TOL) {
        return vec![vec![cx_step(0, 1)]];
    }
    if c.abs() < TOL {
        // CX (Rx(θ) ⊗ Rz(φ)) CX = exp(-iθ/2 XX - iφ/2 ZZ); either coordinate can
        // take either role, and Rx(θ) on the control may be dressed with Z
        // rotations into SX Rz(π - θ) SX.
        let mut out = Vec::new();
        for (on_control, on_target) in [(a, b), (b, a)] {
            for control in [
                rx(-2.0 * on_control),
                sx() * rz(PI - 2.0 * on_control) * sx(),
            ] {
                out.push(vec![
                    cx_step(0, 1),
                    Step::U(0, control),
                    Step::U(1, rz(-2.0 * on_target)),
                    cx_step(0, 1),
                ]);
            }
        }
        return out;
    }
    let half = FRAC_PI_2;
    let mut variants = vec![vec![
        cx_step(1, 0),
        Step::U(0, rz(2.0 * c + half)),
        Step::U(1, ry(2.0 * a + half)),
        cx_step(0, 1),
        Step::U(1, ry(2.0 * b + half)),
        cx_step(1, 0),
    ]];
    // Ry(θ) = Rx(-sπ/2) Rz(sθ) Rx(sπ/2); the X rotations around the middle CX's
    // target cancel through it.
    for sign in [1.0, -1.0] {
        variants.push(vec![
            cx_step(1, 0),
            Step::U(0, rz(2.0 * c + half)),
            Step::U(1, rz(sign * (2.0 * a + half)) * rx(sign * half)),
            cx_step(0, 1),
            Step::U(1, rx(-sign * half) * rz(sign * (2.0 * b + half))),
            cx_step(1, 0),
        ]);
    }
    variants
}

/// A single-qubit Clifford V with V P_src V† = P_dst (axes 0, 1, 2 = X, Y, Z).
fn clifford_map(src: usize, dst: usize) -> M2 {
    use crate::gates::{h, s};
    let (h, s) = (h(), s());
    let candidates = [
        M2::IDENTITY,
        h,
        s,
        s.adjoint(),
        s * h,
        h * s,
        h * s.adjoint(),
        s.adjoint() * h,
    ];
    for v in candidates {
        if (v * pauli(src) * v.adjoint()).close(&pauli(dst), 1e-12) {
            return v;
        }
    }
    unreachable!("no Clifford maps {src} to {dst}")
}

/// exp(i(a XX + b YY + c ZZ)) as one native rotation per non-zero coordinate.
fn ising_template(point: [f64; 3], name: &str, first: usize, second: usize) -> Vec<Step> {
    let mut steps = Vec::new();
    for (axis, coordinate) in point.iter().enumerate() {
        if coordinate.abs() < TOL {
            continue;
        }
        // R_PP(θ) = (V ⊗ W) native(θ) (V ⊗ W)† with V P_first V† = P = W P_second W†.
        let v = clifford_map(first, axis);
        let w = clifford_map(second, axis);
        steps.extend([
            Step::U(0, v.adjoint()),
            Step::U(1, w.adjoint()),
            Step::Gate(name.to_string(), vec![-2.0 * coordinate], [0, 1]),
            Step::U(0, v),
            Step::U(1, w),
        ]);
    }
    steps
}

fn oriented(name: &str, params: &[f64], operands: [usize; 2]) -> M4 {
    let g = matrix_2q(name, params).unwrap_or_else(|| panic!("no matrix for {name}"));
    if operands == [1, 0] { g.swapped() } else { g }
}

/// The unitary of a sequence of steps.
pub fn steps_matrix(steps: &[Step]) -> M4 {
    let mut out = M4::IDENTITY;
    for step in steps {
        let m = match step {
            Step::U(0, m) => kron(m, &M2::IDENTITY),
            Step::U(_, m) => kron(&M2::IDENTITY, m),
            Step::Gate(name, params, q) => oriented(name, params, *q),
        };
        out = m * out;
    }
    out
}

/// Multiply runs of single-qubit unitaries on the same qubit together.
fn merge_locals(steps: Vec<Step>) -> Vec<Step> {
    let mut out = Vec::with_capacity(steps.len());
    let mut pending: [Option<M2>; 2] = [None, None];
    for step in steps {
        match step {
            Step::U(q, m) => pending[q] = Some(m * pending[q].unwrap_or(M2::IDENTITY)),
            Step::Gate(name, params, qubits) => {
                for q in qubits {
                    if let Some(m) = pending[q].take() {
                        out.push(Step::U(q, m));
                    }
                }
                out.push(Step::Gate(name, params, qubits));
            }
        }
    }
    for (q, m) in pending.into_iter().enumerate() {
        if let Some(m) = m {
            out.push(Step::U(q, m));
        }
    }
    out
}

/// Key of a gate for the caches: name, parameter bits and orientation.
type GateKey = (String, Vec<u64>, [usize; 2]);

/// Pauli products `[P, Q]` (0..4 = I, X, Y, Z) a gate maps to Pauli products, with their images.
type PauliImages = Vec<([usize; 2], [usize; 2])>;

/// Two-qubit synthesis with caches of per-gate commutation data.
#[derive(Default)]
pub struct TwoQubitSynth {
    commuting: HashMap<GateKey, [Vec<usize>; 2]>,
    images: HashMap<GateKey, PauliImages>,
}

fn key(name: &str, params: &[f64], operands: [usize; 2]) -> GateKey {
    (
        name.to_string(),
        params.iter().map(|p| p.to_bits()).collect(),
        operands,
    )
}

/// I, X, Y, Z by index 0..4.
fn pauli4(k: usize) -> M2 {
    if k == 0 { M2::IDENTITY } else { pauli(k - 1) }
}

impl TwoQubitSynth {
    pub fn new() -> TwoQubitSynth {
        TwoQubitSynth::default()
    }

    /// `u` as native `basis` gates and single-qubit unitaries, in time order.
    /// With `cost_1q` (the gate count of a single-qubit unitary), the cheapest
    /// template variant is chosen and rotations are moved across the gates.
    pub fn synthesize(
        &mut self,
        u: &M4,
        basis: &TwoQubitBasis,
        cost_1q: Option<&dyn Fn(&M2) -> usize>,
    ) -> Vec<Step> {
        let target = kak(u);
        let templates = match &basis.kind {
            BasisKind::Ising(p, q) => vec![ising_template(target.point, &basis.name, *p, *q)],
            BasisKind::CxLike(_) => cx_templates(target.point),
        };
        let mut candidates: Vec<Vec<Step>> = templates
            .into_iter()
            .filter_map(|t| instantiate(&target, t, basis))
            .collect();
        assert!(
            !candidates.is_empty(),
            "no template reaches the point {:?}",
            target.point
        );
        let Some(cost) = cost_1q else {
            return candidates.swap_remove(0);
        };
        let total = |steps: &[Step]| -> usize {
            steps
                .iter()
                .map(|s| if let Step::U(_, m) = s { cost(m) } else { 1 })
                .sum()
        };
        // Moving rotations is the expensive part: do it for the two variants
        // that start cheapest. The sort is stable, like Python's.
        candidates.sort_by_key(|s| total(s));
        let mut best: Option<(usize, Vec<Step>)> = None;
        for steps in candidates.into_iter().take(2) {
            let dressed = self.dress(steps, cost);
            let t = total(&dressed);
            if best.as_ref().is_none_or(|(b, _)| t < *b) {
                best = Some((t, dressed));
            }
        }
        best.unwrap().1
    }

    fn commuting_axes(
        &mut self,
        name: &str,
        params: &[f64],
        operands: [usize; 2],
    ) -> [Vec<usize>; 2] {
        self.commuting
            .entry(key(name, params, operands))
            .or_insert_with(|| {
                let g = oriented(name, params, operands);
                std::array::from_fn(|q| {
                    (0..3)
                        .filter(|&axis| {
                            let m = if q == 0 {
                                kron(&pauli(axis), &M2::IDENTITY)
                            } else {
                                kron(&M2::IDENTITY, &pauli(axis))
                            };
                            (g * m).close(&(m * g), 1e-12)
                        })
                        .collect()
                })
            })
            .clone()
    }

    fn pauli_images(
        &mut self,
        name: &str,
        params: &[f64],
        operands: [usize; 2],
    ) -> Vec<([usize; 2], [usize; 2])> {
        self.images
            .entry(key(name, params, operands))
            .or_insert_with(|| {
                let g = oriented(name, params, operands);
                let mut out = Vec::new();
                for p in 0..4 {
                    for q in 0..4 {
                        if (p, q) == (0, 0) {
                            continue;
                        }
                        let image = g * kron(&pauli4(p), &pauli4(q)) * g.adjoint();
                        'search: for p2 in 0..4 {
                            for q2 in 0..4 {
                                if (kron(&pauli4(p2), &pauli4(q2)).inner(&image).norm() / 4.0 - 1.0)
                                    .abs()
                                    < 1e-9
                                {
                                    out.push(([p, q], [p2, q2]));
                                    break 'search;
                                }
                            }
                        }
                    }
                }
                out
            })
            .clone()
    }

    /// Shift rotations that commute with a two-qubit gate across it, and Pauli
    /// products across Clifford gates, when that cheapens the neighbours.
    fn dress(&mut self, steps: Vec<Step>, cost: &dyn Fn(&M2) -> usize) -> Vec<Step> {
        let mut work = steps;
        let local_total = |work: &[Step]| -> usize {
            work.iter()
                .map(|s| if let Step::U(_, m) = s { cost(m) } else { 0 })
                .sum()
        };
        for _ in 0..3 {
            let total = local_total(&work);
            let mut i = 0;
            while i < work.len() {
                let (name, params, operands) = match &work[i] {
                    Step::U(..) => {
                        i += 1;
                        continue;
                    }
                    Step::Gate(name, params, operands) => (name.clone(), params.clone(), *operands),
                };
                let axes = self.commuting_axes(&name, &params, operands);
                for q in operands {
                    for &axis in &axes[q] {
                        let before = match neighbour(&work, i, q, -1) {
                            Some(b) => b,
                            None => {
                                work.insert(i, Step::U(q, M2::IDENTITY));
                                i += 1;
                                i - 1
                            }
                        };
                        let after = match neighbour(&work, i, q, 1) {
                            Some(a) => a,
                            None => {
                                work.insert(i + 1, Step::U(q, M2::IDENTITY));
                                i + 1
                            }
                        };
                        let (b, a) = (local(&work[before]), local(&work[after]));
                        let mut best_cost = cost(&b) + cost(&a);
                        let mut best_angle = None;
                        for angle in candidate_angles(axis, &b, &a) {
                            let r = rotation(axis, angle);
                            let c = cost(&(r.adjoint() * b)) + cost(&(a * r));
                            if c < best_cost {
                                (best_cost, best_angle) = (c, Some(angle));
                            }
                        }
                        if let Some(angle) = best_angle {
                            let r = rotation(axis, angle);
                            set_local(&mut work[before], r.adjoint() * b);
                            set_local(&mut work[after], a * r);
                        }
                    }
                }
                i = self.move_paulis(&mut work, i, cost);
                i += 1;
            }
            if local_total(&work) >= total {
                break;
            }
        }
        work.retain(|s| !matches!(s, Step::U(_, m) if m.close(&M2::IDENTITY, 1e-8)));
        work
    }

    /// Move a Pauli product P⊗Q across Clifford gate `i` if that cheapens its
    /// neighbours: G (P⊗Q) = (P'⊗Q') G. Returns the (possibly shifted) index of the gate.
    fn move_paulis(
        &mut self,
        work: &mut Vec<Step>,
        mut i: usize,
        cost: &dyn Fn(&M2) -> usize,
    ) -> usize {
        let Step::Gate(name, params, operands) = &work[i] else {
            return i;
        };
        let images = self.pauli_images(&name.clone(), &params.clone(), *operands);
        if images.is_empty() {
            return i;
        }
        // Make sure both qubits have a single-qubit step on each side, then find them.
        for q in 0..2 {
            if neighbour(work, i, q, -1).is_none() {
                work.insert(i, Step::U(q, M2::IDENTITY));
                i += 1;
            }
            if neighbour(work, i, q, 1).is_none() {
                work.insert(i + 1, Step::U(q, M2::IDENTITY));
            }
        }
        let slots = [
            neighbour(work, i, 0, -1).unwrap(),
            neighbour(work, i, 1, -1).unwrap(),
            neighbour(work, i, 0, 1).unwrap(),
            neighbour(work, i, 1, 1).unwrap(),
        ];
        let [b0, b1, a0, a1] = slots.map(|j| local(&work[j]));
        let mut best_cost = cost(&b0) + cost(&b1) + cost(&a0) + cost(&a1);
        let mut best = None;
        for ([p, q], [p2, q2]) in images {
            let moved = [
                pauli4(p) * b0,
                pauli4(q) * b1,
                a0 * pauli4(p2),
                a1 * pauli4(q2),
            ];
            let c: usize = moved.iter().map(cost).sum();
            if c < best_cost {
                (best_cost, best) = (c, Some(moved));
            }
        }
        if let Some(moved) = best {
            for (j, m) in slots.into_iter().zip(moved) {
                set_local(&mut work[j], m);
            }
        }
        i
    }
}

fn local(step: &Step) -> M2 {
    match step {
        Step::U(_, m) => *m,
        Step::Gate(..) => unreachable!(),
    }
}

fn set_local(step: &mut Step, value: M2) {
    if let Step::U(_, m) = step {
        *m = value;
    }
}

/// The index of the single-qubit step on `q` next to step `i`, if no gate on `q` intervenes.
fn neighbour(steps: &[Step], i: usize, q: usize, direction: isize) -> Option<usize> {
    let mut j = i as isize + direction;
    while 0 <= j && (j as usize) < steps.len() {
        match &steps[j as usize] {
            Step::U(on, _) if *on == q => return Some(j as usize),
            Step::Gate(_, _, operands) if operands.contains(&q) => return None,
            _ => {}
        }
        j += direction;
    }
    None
}

/// Angles α where `f(α) = |m(α)_00|²` is 1, 1/2 or 0, i.e. where the Euler
/// angle θ of `m(α)` is 0, π/2 or π, which many bases implement with fewer gates.
fn special_angles(f: impl Fn(f64) -> f64) -> Vec<f64> {
    let (f0, f1, f2) = (f(0.0), f(FRAC_PI_2), f(PI));
    let (a, b, c) = ((f0 + f2) / 2.0, (f0 - f2) / 2.0, f1 - (f0 + f2) / 2.0);
    let r = b.hypot(c);
    if r < 1e-12 {
        return vec![];
    }
    let delta = c.atan2(b);
    let mut out = Vec::new();
    for target in [1.0, 0.5, 0.0] {
        let k = (target - a) / r;
        if k.abs() <= 1.0 + 1e-12 {
            let x = k.clamp(-1.0, 1.0).acos();
            out.push(delta + x);
            out.push(delta - x);
        }
    }
    out
}

/// Angles worth trying for a rotation about `axis` moved from `before` to `after`.
fn candidate_angles(axis: usize, before: &M2, after: &M2) -> Vec<f64> {
    let mut out = vec![0.0];
    out.extend(special_angles(|t| {
        (rotation(axis, t).adjoint() * *before).0[0][0].norm_sqr()
    }));
    out.extend(special_angles(|t| {
        (*after * rotation(axis, t)).0[0][0].norm_sqr()
    }));
    if axis == 2 {
        // Z rotations change only the outer Euler angles: cancel one of them,
        // or both together when θ = 0.
        let (_, phi, lam_b) = euler_zyz(before);
        let (_, phi_a, lam) = euler_zyz(after);
        for base in [phi, -lam, phi + lam_b, -(phi_a + lam)] {
            for extra in [0.0, PI, FRAC_PI_2, -FRAC_PI_2] {
                out.push(base + extra);
            }
        }
    }
    // Dedupe modulo 2π, keeping the first position and the last value, like a
    // Python dict.
    let mut keys: Vec<i64> = Vec::new();
    let mut unique: Vec<f64> = Vec::new();
    for angle in out {
        let k = (remainder(angle, TAU) * 1e9).round() as i64;
        match keys.iter().position(|&x| x == k) {
            Some(pos) => unique[pos] = angle,
            None => {
                keys.push(k);
                unique.push(angle);
            }
        }
    }
    unique
}

/// `target` as `template` between single-qubit gates, or None if the points differ.
fn instantiate(target: &Kak, template: Vec<Step>, basis: &TwoQubitBasis) -> Option<Vec<Step>> {
    let reference = kak(&steps_matrix(&template));
    if !close3(reference.point, target.point, 1e-7) {
        return None;
    }
    // u ∝ (A C†) T (D† B) on each qubit, where T is the template.
    let mut steps = vec![
        Step::U(0, reference.b1.adjoint() * target.b1),
        Step::U(1, reference.b0.adjoint() * target.b0),
    ];
    for step in template {
        match (&step, &basis.kind) {
            (Step::Gate(name, _, q), BasisKind::CxLike(locals))
                if name == "cx" && basis.name != "cx" =>
            {
                let [after1, after0, before1, before0] = **locals;
                let [control, target_qubit] = *q;
                steps.extend([
                    Step::U(control, before1),
                    Step::U(target_qubit, before0),
                    Step::Gate(basis.name.clone(), vec![], *q),
                    Step::U(control, after1),
                    Step::U(target_qubit, after0),
                ]);
            }
            _ => steps.push(step),
        }
    }
    steps.push(Step::U(0, target.a1 * reference.a1.adjoint()));
    steps.push(Step::U(1, target.a0 * reference.a0.adjoint()));
    Some(merge_locals(steps))
}

/// Whether `steps` implement `u` up to a global phase.
pub fn implements(steps: &[Step], u: &M4) -> bool {
    equal_up_to_phase4(&steps_matrix(steps), u, 1e-8)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gates::{matrix_2q, signature};
    use crate::one_qubit::{bases_for, count_1q};
    use crate::testing::Rng;

    fn standard_2q() -> Vec<M4> {
        let names = [
            "cx", "cy", "cz", "ch", "csx", "swap", "iswap", "dcx", "ecr", "cp", "crx", "cry",
            "crz", "rxx", "ryy", "rzz", "rzx", "cu3", "cu",
        ];
        names
            .iter()
            .map(|n| {
                let count = signature(n).unwrap().1;
                let params: Vec<f64> = (0..count).map(|i| 0.37 + 0.41 * i as f64).collect();
                matrix_2q(n, &params).unwrap()
            })
            .collect()
    }

    #[test]
    fn kak_reconstructs_and_is_canonical() {
        let mut rng = Rng::new(7);
        for _ in 0..2000 {
            let u = rng.unitary4();
            let k = kak(&u);
            let [a, b, c] = k.point;
            assert!(
                FRAC_PI_4 + 1e-9 >= a && a >= b - 1e-9 && b >= c.abs() - 1e-9,
                "{:?}",
                k.point
            );
            assert!(k.matrix().close(&u, 1e-9));
        }
        for u in standard_2q() {
            assert!(kak(&u).matrix().close(&u, 1e-9));
        }
    }

    #[test]
    fn weyl_points() {
        let q = FRAC_PI_4;
        for (name, point) in [
            ("cx", [q, 0.0, 0.0]),
            ("cz", [q, 0.0, 0.0]),
            ("ecr", [q, 0.0, 0.0]),
            ("iswap", [q, q, 0.0]),
            ("swap", [q, q, q]),
        ] {
            let p = kak(&matrix_2q(name, &[]).unwrap()).point;
            assert!(
                (0..3).all(|i| (p[i] - point[i]).abs() < 1e-9),
                "{name} {p:?}"
            );
        }
    }

    #[test]
    fn synthesis_is_exact_with_and_without_dressing() {
        let bases_1q = bases_for(&["rz", "sx", "x"].iter().map(|s| s.to_string()).collect());
        let cost = |m: &M2| count_1q(m, &bases_1q);
        let mut rng = Rng::new(5);
        let mut samples: Vec<M4> = (0..100).map(|_| rng.unitary4()).collect();
        samples.extend(standard_2q());
        let mut synth = TwoQubitSynth::new();
        for native in ["cx", "cz", "cy", "ch", "ecr", "rzz", "rxx", "ryy", "rzx"] {
            let basis = two_qubit_basis(native).unwrap();
            for u in &samples {
                for steps in [
                    synth.synthesize(u, &basis, None),
                    synth.synthesize(u, &basis, Some(&cost)),
                ] {
                    assert!(implements(&steps, u), "{native}");
                    assert!(steps.iter().all(|s| matches!(s, Step::U(..))
                        || matches!(s, Step::Gate(n, ..) if n == native)));
                }
            }
        }
    }

    #[test]
    fn cx_counts_are_optimal() {
        let basis = two_qubit_basis("cx").unwrap();
        let mut synth = TwoQubitSynth::new();
        for (name, params, count) in [
            ("cz", vec![], 1),
            ("ch", vec![], 1),
            ("crz", vec![0.37], 2),
            ("iswap", vec![], 2),
            ("swap", vec![], 3),
        ] {
            let steps = synth.synthesize(&matrix_2q(name, &params).unwrap(), &basis, None);
            assert_eq!(
                steps.iter().filter(|s| matches!(s, Step::Gate(..))).count(),
                count,
                "{name}"
            );
        }
        assert!(two_qubit_basis("iswap").is_none());
    }
}
