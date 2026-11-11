//! Optimal single-qubit synthesis into a device's native gates.
//!
//! Every single-qubit unitary is, up to global phase, `Rz(φ) Ry(θ) Rz(λ)` (a
//! ZYZ Euler decomposition). Each native basis turns those angles into the
//! fewest gates it allows, dropping rotations by multiples of 2π. Matches
//! `qevedo/compiler/synthesis/one_qubit.py`.

use crate::gates::h;
use crate::ir::{Inline, Name};
use crate::linalg::{C, M2};
use std::collections::HashSet;
use std::f64::consts::{FRAC_PI_2, PI, TAU};

pub const TOL: f64 = 1e-10;

/// One gate of a synthesized sequence.
#[derive(Clone, Debug, PartialEq)]
pub struct Op {
    pub name: Name,
    pub params: Inline<f64, 3>,
}

impl Op {
    fn new(name: &str, params: impl Into<Inline<f64, 3>>) -> Op {
        Op {
            name: Name::from(name),
            params: params.into(),
        }
    }
}

/// `x` modulo 2π, in [-π, π], rounding half to even like Python's `math.remainder`.
pub fn remainder(x: f64, period: f64) -> f64 {
    x - (x / period).round_ties_even() * period
}

/// `angle` in (-π, π].
pub fn wrap(angle: f64) -> f64 {
    let wrapped = remainder(angle, TAU);
    if (wrapped + PI).abs() <= TOL {
        PI
    } else {
        wrapped
    }
}

pub fn is_zero(angle: f64) -> bool {
    wrap(angle).abs() < TOL
}

/// ZYZ angles `(θ, φ, λ)` with `u ≅ Rz(φ) Ry(θ) Rz(λ)` up to global phase, θ in [0, π].
pub fn euler_zyz(u: &M2) -> (f64, f64, f64) {
    let su = u.scale(u.det().sqrt().inv());
    // su = [[e^{-i(φ+λ)/2} cos(θ/2), …], [e^{i(φ-λ)/2} sin(θ/2), e^{i(φ+λ)/2} cos(θ/2)]]
    let (a, b) = (su.0[1][1], su.0[1][0]);
    let theta = 2.0 * b.norm().atan2(a.norm());
    let plus = if a.norm() > TOL { 2.0 * a.arg() } else { 0.0 };
    let minus = if b.norm() > TOL { 2.0 * b.arg() } else { 0.0 };
    (theta, (plus + minus) / 2.0, (plus - minus) / 2.0)
}

/// A family of native single-qubit gates. `z` names the gate used as Rz (it may
/// be `p` or `u1`, which differ from `rz` only by a global phase).
#[derive(Clone, Debug)]
pub enum Basis1q {
    U3(String),
    Zsx { z: String, x: bool },
    Zxz { z: String },
    Zyz { z: String },
    Xyx,
}

fn rot(name: &str, angle: f64, out: &mut Vec<Op>) {
    if !is_zero(angle) {
        out.push(Op::new(name, [wrap(angle)]));
    }
}

fn nonzero(angles: &[f64]) -> usize {
    angles.iter().filter(|a| !is_zero(**a)).count()
}

impl Basis1q {
    /// The shortest sequence of this basis (in time order) equal to `u` up to phase.
    pub fn emit(&self, u: &M2) -> Vec<Op> {
        let mut out = Vec::with_capacity(5);
        match self {
            Basis1q::U3(name) => {
                let (theta, phi, lam) = euler_zyz(u);
                if !(is_zero(theta) && is_zero(phi + lam)) {
                    out.push(Op::new(name, [wrap(theta), wrap(phi), wrap(lam)]));
                }
            }
            Basis1q::Zsx { z, x } => {
                let (theta, phi, lam) = euler_zyz(u);
                if is_zero(theta) {
                    rot(z, phi + lam, &mut out);
                } else if (theta - FRAC_PI_2).abs() < TOL {
                    // U3(π/2, φ, λ) ≅ Rz(φ + π/2) SX Rz(λ - π/2)
                    rot(z, lam - FRAC_PI_2, &mut out);
                    out.push(Op::new("sx", []));
                    rot(z, phi + FRAC_PI_2, &mut out);
                } else if *x && (theta - PI).abs() < TOL {
                    // U3(π, φ, λ) ≅ Rz(φ + π) X Rz(λ) = Rz(φ - λ + π) X
                    out.push(Op::new("x", []));
                    rot(z, phi - lam + PI, &mut out);
                } else {
                    // U3(θ, φ, λ) ≅ Rz(φ + π) SX Rz(θ + π) SX Rz(λ)
                    rot(z, lam, &mut out);
                    out.push(Op::new("sx", []));
                    rot(z, theta + PI, &mut out);
                    out.push(Op::new("sx", []));
                    rot(z, phi + PI, &mut out);
                }
            }
            Basis1q::Zxz { z } => {
                // Ry(θ) = Rz(π/2) Rx(θ) Rz(-π/2)
                let (theta, phi, lam) = euler_zyz(u);
                if is_zero(theta) {
                    rot(z, phi + lam, &mut out);
                } else {
                    rot(z, lam - FRAC_PI_2, &mut out);
                    out.push(Op::new("rx", [wrap(theta)]));
                    rot(z, phi + FRAC_PI_2, &mut out);
                }
            }
            Basis1q::Zyz { z } => {
                let (theta, phi, lam) = euler_zyz(u);
                if is_zero(theta) {
                    rot(z, phi + lam, &mut out);
                } else {
                    rot(z, lam, &mut out);
                    out.push(Op::new("ry", [wrap(theta)]));
                    rot(z, phi, &mut out);
                }
            }
            Basis1q::Xyx => {
                // H Rz(a) H = Rx(a) and H Ry(b) H = Ry(-b): the ZYZ angles of H U H.
                let (theta, phi, lam) = euler_zyz(&(h() * *u * h()));
                if is_zero(theta) {
                    rot("rx", phi + lam, &mut out);
                } else {
                    rot("rx", lam, &mut out);
                    out.push(Op::new("ry", [wrap(-theta)]));
                    rot("rx", phi, &mut out);
                }
            }
        }
        out
    }

    /// `self.emit(u).len()`, without building the sequence, from what
    /// [`Shape`] reads off `u`'s entries; from the ZYZ angles
    /// ([`Basis1q::count_from_angles`]) only where the shape cannot tell.
    pub fn count(&self, u: &M2, shape: &Shape) -> usize {
        let fallback = || self.count_from_angles(u, euler_zyz(u));
        match self {
            Basis1q::U3(_) => usize::from(!(shape.theta_zero && is_one(shape.plus))),
            Basis1q::Zsx { x, .. } => {
                if shape.theta_zero {
                    nonzero_turns(&[shape.plus])
                } else if shape.theta_half {
                    // λ - π/2 and φ + π/2
                    match (shape.lambda, shape.phi) {
                        (Some(l), Some(p)) => 1 + nonzero_turns(&[l * -I, p * I]),
                        _ => fallback(),
                    }
                } else if *x && shape.theta_pi {
                    // φ - λ + π
                    1 + nonzero_turns(&[-shape.minus])
                } else {
                    // λ, θ + π (zero exactly when θ = π) and φ + π
                    match (shape.lambda, shape.phi) {
                        (Some(l), Some(p)) => {
                            2 + nonzero_turns(&[l, -p]) + usize::from(!shape.theta_pi)
                        }
                        _ => fallback(),
                    }
                }
            }
            Basis1q::Zxz { .. } | Basis1q::Zyz { .. } => {
                if shape.theta_zero {
                    return nonzero_turns(&[shape.plus]);
                }
                let (Some(l), Some(p)) = (shape.lambda, shape.phi) else {
                    return fallback();
                };
                match self {
                    // λ - π/2 and φ + π/2
                    Basis1q::Zxz { .. } => 1 + nonzero_turns(&[l * -I, p * I]),
                    _ => 1 + nonzero_turns(&[l, p]),
                }
            }
            Basis1q::Xyx => self.emit(u).len(),
        }
    }

    /// `self.emit(u).len()`, given the ZYZ angles of `u`, without building the sequence.
    pub fn count_from_angles(&self, u: &M2, angles: (f64, f64, f64)) -> usize {
        let (theta, phi, lam) = angles;
        match self {
            Basis1q::U3(_) => usize::from(!(is_zero(theta) && is_zero(phi + lam))),
            Basis1q::Zsx { x, .. } => {
                if is_zero(theta) {
                    nonzero(&[phi + lam])
                } else if (theta - FRAC_PI_2).abs() < TOL {
                    1 + nonzero(&[lam - FRAC_PI_2, phi + FRAC_PI_2])
                } else if *x && (theta - PI).abs() < TOL {
                    1 + nonzero(&[phi - lam + PI])
                } else {
                    2 + nonzero(&[lam, theta + PI, phi + PI])
                }
            }
            Basis1q::Zxz { .. } => {
                if is_zero(theta) {
                    nonzero(&[phi + lam])
                } else {
                    1 + nonzero(&[lam - FRAC_PI_2, phi + FRAC_PI_2])
                }
            }
            Basis1q::Zyz { .. } => {
                if is_zero(theta) {
                    nonzero(&[phi + lam])
                } else {
                    1 + nonzero(&[lam, phi])
                }
            }
            Basis1q::Xyx => self.emit(u).len(),
        }
    }
}

/// Gates that act as Rz up to a global phase, in order of preference.
const Z_GATES: [&str; 4] = ["rz", "p", "u1", "phase"];

/// The single-qubit bases a native gate set supports.
pub fn bases_for(native: &HashSet<String>) -> Vec<Basis1q> {
    let z = Z_GATES
        .iter()
        .find(|g| native.contains(**g))
        .map(|g| g.to_string());
    let mut out = Vec::new();
    if let Some(name) = ["u3", "u"].iter().find(|g| native.contains(**g)) {
        out.push(Basis1q::U3(name.to_string()));
    }
    if let Some(z) = &z {
        if native.contains("sx") {
            out.push(Basis1q::Zsx {
                z: z.clone(),
                x: native.contains("x"),
            });
        }
        if native.contains("rx") {
            out.push(Basis1q::Zxz { z: z.clone() });
        }
        if native.contains("ry") {
            out.push(Basis1q::Zyz { z: z.clone() });
        }
    }
    if native.contains("rx") && native.contains("ry") {
        out.push(Basis1q::Xyx);
    }
    out
}

/// The shortest sequence over `bases` (in time order) equal to `u` up to phase.
pub fn synthesize_1q(u: &M2, bases: &[Basis1q]) -> Vec<Op> {
    let mut best: Option<Vec<Op>> = None;
    for basis in bases {
        let ops = basis.emit(u);
        if best.as_ref().is_none_or(|b| ops.len() < b.len()) {
            best = Some(ops);
        }
    }
    best.expect("the native gate set has no universal single-qubit basis")
}

/// `synthesize_1q(u, bases).len()`, without building the sequence.
pub fn count_1q(u: &M2, bases: &[Basis1q]) -> usize {
    let shape = Shape::of(u);
    bases.iter().map(|b| b.count(u, &shape)).min().unwrap_or(0)
}

/// What [`Basis1q::count`] needs of the ZYZ angles of a unitary, read off
/// its entries with products and square roots instead of the inverse
/// trigonometric functions [`euler_zyz`] calls.
///
/// Counting is what two-qubit synthesis does most: it prices every way of
/// moving a rotation across a gate by the single-qubit gates it leaves, and
/// on the 32-qubit Quantum Volume circuit it counts some 650 000 unitaries,
/// where `atan2`, `hypot` and the complex square root of the determinant were
/// a quarter of the compilation. A count only asks whether some angle is a
/// multiple of 2π, and the angles are those of products of the entries: with
/// `s = det(u)^(-1/2)`, `a = s·u₁₁ = e^{i(φ+λ)/2} cos(θ/2)` and
/// `b = s·u₁₀ = e^{i(φ-λ)/2} sin(θ/2)`,
///
/// - `e^{i(φ+λ)} ∝ a² ∝ u₁₁² · conj(det u)` and `e^{i(φ-λ)} ∝ u₁₀² · conj(det u)`,
/// - `e^{iφ} ∝ a·b ∝ u₁₁ u₁₀ · conj(det u)` and `e^{iλ} ∝ a·conj(b) ∝ u₁₁ · conj(u₁₀)`,
///
/// up to positive factors, which do not change an angle. The square root's
/// sign never enters (s appears squared, or times its conjugate). θ is read
/// from |a| and |b| as `euler_zyz` does, compared through its tangent or
/// cosine rather than computed.
///
/// [`euler_zyz`] drops the phase of an entry no larger than [`TOL`] (its
/// angle is then 0). When exactly one of `a` and `b` is that small, φ and λ
/// alone are half that angle, whose direction needs `s` itself, and `phi`
/// and `lambda` are None: the count then falls back
/// to the angles. The sums the bases need there (φ + λ when θ = 0, φ - λ when
/// θ = π) do not.
pub struct Shape {
    /// θ = 0, π/2 and π, within [`TOL`].
    theta_zero: bool,
    theta_half: bool,
    theta_pi: bool,
    /// Proportional to e^{i(φ+λ)} and e^{i(φ-λ)} (1 where `euler_zyz` drops the angle).
    plus: C,
    minus: C,
    /// Proportional to e^{iφ} and e^{iλ}, unless exactly one of `a` and `b`
    /// loses its phase.
    phi: Option<C>,
    lambda: Option<C>,
}

impl Shape {
    pub fn of(u: &M2) -> Shape {
        let (a, b) = (u.0[1][1], u.0[1][0]);
        let det = u.det();
        // |s| = |det u|^(-1/2) is 1 for a unitary, up to rounding.
        let scale = det.norm_sqr().sqrt().sqrt();
        let (na2, nb2) = (a.norm_sqr(), b.norm_sqr());
        let (na, nb) = (na2.sqrt(), nb2.sqrt());
        let (keep_a, keep_b) = (na > TOL * scale, nb > TOL * scale);
        let conj_det = det.conj();
        // θ/2 = atan2(|b|, |a|) is below TOL/2 when |b| < tan(TOL/2)·|a|,
        // and θ is within TOL of π/2 when |cos θ| < sin(TOL).
        Shape {
            theta_zero: nb < TAN_HALF_TOL * na,
            theta_half: (na2 - nb2).abs() < SIN_TOL * (na2 + nb2),
            theta_pi: na < TAN_HALF_TOL * nb,
            plus: if keep_a { a * a * conj_det } else { ONE_C },
            minus: if keep_b { b * b * conj_det } else { ONE_C },
            phi: match (keep_a, keep_b) {
                (true, true) => Some(a * b * conj_det),
                (false, false) => Some(ONE_C),
                _ => None,
            },
            lambda: match (keep_a, keep_b) {
                (true, true) => Some(a * b.conj()),
                (false, false) => Some(ONE_C),
                _ => None,
            },
        }
    }
}

const ONE_C: C = C::new(1.0, 0.0);
/// tan(TOL/2) and sin(TOL), which differ from TOL/2 and TOL by parts in
/// 10²¹, far below a double's precision.
const TAN_HALF_TOL: f64 = TOL / 2.0;
const SIN_TOL: f64 = TOL;
const I: C = C::new(0.0, 1.0);

/// Whether the angle of `w` (any non-zero multiple of `e^{iα}` by a positive
/// number) is a multiple of 2π within [`TOL`], as [`is_zero`] says of α:
/// |α| < TOL modulo 2π exactly when `w` points right and |sin α| < sin(TOL).
fn is_one(w: C) -> bool {
    w.re > 0.0 && w.im.abs() < SIN_TOL * w.norm_sqr().sqrt()
}

/// How many of the angles of `turns` (see [`is_one`]) are not multiples of 2π.
fn nonzero_turns(turns: &[C]) -> usize {
    turns.iter().filter(|w| !is_one(**w)).count()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gates::matrix_1q;
    use crate::testing::{RandomUnitary, Rng, equal_up_to_phase2};

    fn ops_matrix(ops: &[Op]) -> M2 {
        ops.iter().fold(M2::IDENTITY, |m, op| {
            matrix_1q(&op.name, &op.params).unwrap() * m
        })
    }

    fn native(names: &[&str]) -> HashSet<String> {
        names.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn every_basis_is_exact_and_counts_match() {
        let sets: [&[&str]; 7] = [
            &["rz", "sx", "x"],
            &["rz", "sx"],
            &["rz", "rx"],
            &["rz", "ry"],
            &["rx", "ry"],
            &["u3"],
            &["p", "sx"],
        ];
        let mut rng = Rng::new(3);
        for set in sets {
            let bases = bases_for(&native(set));
            let mut samples: Vec<M2> = (0..300).map(|_| rng.unitary2()).collect();
            for name in ["id", "x", "y", "z", "h", "s", "t", "sx", "sxdg"] {
                samples.push(matrix_1q(name, &[]).unwrap());
            }
            for u in samples {
                let ops = synthesize_1q(&u, &bases);
                assert!(
                    equal_up_to_phase2(&ops_matrix(&ops), &u, 1e-9),
                    "{set:?} {ops:?}"
                );
                assert_eq!(count_1q(&u, &bases), ops.len());
                assert!(ops.iter().all(|op| set.contains(&op.name.as_str())));
            }
        }
    }

    #[test]
    fn counting_from_the_entries_agrees_with_the_angles() {
        // Unitaries at the special angles a count looks for, as two-qubit
        // synthesis makes them (θ of 0, π/2 or π, outer angles that are
        // multiples of π/2, entries that vanish), times random phases and
        // nudged to either side of the tolerance (not onto it, where
        // rounding decides either way), besides random ones.
        let mut rng = Rng::new(11);
        let quarter = |k: usize| k as f64 * FRAC_PI_2;
        let mut samples: Vec<M2> = (0..300).map(|_| rng.unitary2()).collect();
        for t in 0..5 {
            for p in 0..4 {
                for l in 0..4 {
                    for nudge in [0.0, 1e-13, 2e-11, 3e-10, 1e-9] {
                        let theta = t as f64 * PI / 4.0 + nudge;
                        let (phi, lam) = (quarter(p) - nudge, quarter(l) + 2.0 * nudge);
                        let u = matrix_1q("u3", &[theta, phi, lam]).unwrap();
                        let phase = C::from_polar(1.0, rng.uniform() * TAU);
                        samples.push(u.scale(phase));
                        samples.push(u);
                    }
                }
            }
        }
        let sets: [&[&str]; 6] = [
            &["rz", "sx", "x"],
            &["rz", "sx"],
            &["rz", "rx"],
            &["rz", "ry"],
            &["rx", "ry"],
            &["u3"],
        ];
        for set in sets {
            for basis in bases_for(&native(set)) {
                for u in &samples {
                    assert_eq!(
                        basis.count(u, &Shape::of(u)),
                        basis.count_from_angles(u, euler_zyz(u)),
                        "{basis:?} {u:?}"
                    );
                    assert_eq!(basis.count(u, &Shape::of(u)), basis.emit(u).len());
                }
            }
        }
    }

    #[test]
    fn special_cases() {
        let zsx = bases_for(&native(&["rz", "sx", "x"]));
        let names = |u: M2| {
            synthesize_1q(&u, &zsx)
                .into_iter()
                .map(|o| o.name)
                .collect::<Vec<_>>()
        };
        assert!(names(M2::IDENTITY).is_empty());
        assert_eq!(names(h()), ["rz", "sx", "rz"]);
        assert_eq!(names(crate::gates::x()), ["x"]);
        assert_eq!(names(crate::gates::y()), ["x", "rz"]);
        assert_eq!(names(crate::gates::sx()), ["sx"]);
    }
}
