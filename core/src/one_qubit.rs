//! Optimal single-qubit synthesis into a device's native gates.
//!
//! Every single-qubit unitary is, up to global phase, `Rz(φ) Ry(θ) Rz(λ)` (a
//! ZYZ Euler decomposition). Each native basis turns those angles into the
//! fewest gates it allows, dropping rotations by multiples of 2π. Matches
//! `qevedo/compiler/synthesis/one_qubit.py`.

use crate::gates::h;
use crate::linalg::M2;
use std::collections::HashSet;
use std::f64::consts::{FRAC_PI_2, PI, TAU};

pub const TOL: f64 = 1e-10;

/// One gate of a synthesized sequence.
#[derive(Clone, Debug, PartialEq)]
pub struct Op {
    pub name: String,
    pub params: Vec<f64>,
}

impl Op {
    fn new(name: &str, params: Vec<f64>) -> Op {
        Op {
            name: name.to_string(),
            params,
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
        out.push(Op::new(name, vec![wrap(angle)]));
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
                    out.push(Op::new(name, vec![wrap(theta), wrap(phi), wrap(lam)]));
                }
            }
            Basis1q::Zsx { z, x } => {
                let (theta, phi, lam) = euler_zyz(u);
                if is_zero(theta) {
                    rot(z, phi + lam, &mut out);
                } else if (theta - FRAC_PI_2).abs() < TOL {
                    // U3(π/2, φ, λ) ≅ Rz(φ + π/2) SX Rz(λ - π/2)
                    rot(z, lam - FRAC_PI_2, &mut out);
                    out.push(Op::new("sx", vec![]));
                    rot(z, phi + FRAC_PI_2, &mut out);
                } else if *x && (theta - PI).abs() < TOL {
                    // U3(π, φ, λ) ≅ Rz(φ + π) X Rz(λ) = Rz(φ - λ + π) X
                    out.push(Op::new("x", vec![]));
                    rot(z, phi - lam + PI, &mut out);
                } else {
                    // U3(θ, φ, λ) ≅ Rz(φ + π) SX Rz(θ + π) SX Rz(λ)
                    rot(z, lam, &mut out);
                    out.push(Op::new("sx", vec![]));
                    rot(z, theta + PI, &mut out);
                    out.push(Op::new("sx", vec![]));
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
                    out.push(Op::new("rx", vec![wrap(theta)]));
                    rot(z, phi + FRAC_PI_2, &mut out);
                }
            }
            Basis1q::Zyz { z } => {
                let (theta, phi, lam) = euler_zyz(u);
                if is_zero(theta) {
                    rot(z, phi + lam, &mut out);
                } else {
                    rot(z, lam, &mut out);
                    out.push(Op::new("ry", vec![wrap(theta)]));
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
                    out.push(Op::new("ry", vec![wrap(-theta)]));
                    rot("rx", phi, &mut out);
                }
            }
        }
        out
    }

    /// `self.emit(u).len()`, given the ZYZ angles of `u`, without building the sequence.
    pub fn count(&self, u: &M2, angles: (f64, f64, f64)) -> usize {
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
    let angles = euler_zyz(u);
    bases.iter().map(|b| b.count(u, angles)).min().unwrap_or(0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gates::matrix_1q;
    use crate::testing::{Rng, equal_up_to_phase2};

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
