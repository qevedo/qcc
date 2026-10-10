//! Unitary matrices of the standard OpenQASM gates (`qelib1.inc`,
//! `stdgates.inc`, and Qiskit for `ecr`, `iswap` and `rzx`), with their global
//! phases. Matches `qevedo/compiler/synthesis/gates.py`.

use crate::linalg::{C, I, M2, M4, ONE, ZERO, kron};
use std::f64::consts::{FRAC_1_SQRT_2, FRAC_PI_2, FRAC_PI_4};

pub fn x() -> M2 {
    M2::new(ZERO, ONE, ONE, ZERO)
}
pub fn y() -> M2 {
    M2::new(ZERO, -I, I, ZERO)
}
pub fn z() -> M2 {
    M2::new(ONE, ZERO, ZERO, -ONE)
}
pub fn h() -> M2 {
    let r = C::new(FRAC_1_SQRT_2, 0.0);
    M2::new(r, r, r, -r)
}
pub fn s() -> M2 {
    M2::new(ONE, ZERO, ZERO, I)
}
pub fn sx() -> M2 {
    let (p, m) = (C::new(0.5, 0.5), C::new(0.5, -0.5));
    M2::new(p, m, m, p)
}
pub fn t() -> M2 {
    M2::new(ONE, ZERO, ZERO, C::from_polar(1.0, FRAC_PI_4))
}

/// The Paulis X, Y, Z by axis 0, 1, 2.
pub fn pauli(axis: usize) -> M2 {
    match axis {
        0 => x(),
        1 => y(),
        _ => z(),
    }
}

pub fn rx(theta: f64) -> M2 {
    let (c, s) = ((theta / 2.0).cos(), (theta / 2.0).sin());
    M2::new(
        C::new(c, 0.0),
        C::new(0.0, -s),
        C::new(0.0, -s),
        C::new(c, 0.0),
    )
}
pub fn ry(theta: f64) -> M2 {
    let (c, s) = ((theta / 2.0).cos(), (theta / 2.0).sin());
    M2::new(
        C::new(c, 0.0),
        C::new(-s, 0.0),
        C::new(s, 0.0),
        C::new(c, 0.0),
    )
}
pub fn rz(theta: f64) -> M2 {
    M2::new(
        C::from_polar(1.0, -theta / 2.0),
        ZERO,
        ZERO,
        C::from_polar(1.0, theta / 2.0),
    )
}
pub fn phase(lam: f64) -> M2 {
    M2::new(ONE, ZERO, ZERO, C::from_polar(1.0, lam))
}
pub fn u3(theta: f64, phi: f64, lam: f64) -> M2 {
    let (c, s) = ((theta / 2.0).cos(), (theta / 2.0).sin());
    M2::new(
        C::new(c, 0.0),
        -C::from_polar(s, lam),
        C::from_polar(s, phi),
        C::from_polar(c, phi + lam),
    )
}

/// exp(-iθ/2 P) for the Pauli on `axis`.
pub fn rotation(axis: usize, theta: f64) -> M2 {
    match axis {
        0 => rx(theta),
        1 => ry(theta),
        _ => rz(theta),
    }
}

/// `u` controlled by the first qubit.
pub fn controlled(u: &M2) -> M4 {
    let mut out = M4::IDENTITY;
    out.0[2][2] = u.0[0][0];
    out.0[2][3] = u.0[0][1];
    out.0[3][2] = u.0[1][0];
    out.0[3][3] = u.0[1][1];
    out
}

/// exp(-iθ/2 P⊗Q).
pub fn pauli_rotation(p: usize, q: usize, theta: f64) -> M4 {
    let pq = kron(&pauli(p), &pauli(q));
    let (c, s) = ((theta / 2.0).cos(), (theta / 2.0).sin());
    let mut out = pq.scale(C::new(0.0, -s));
    for i in 0..4 {
        out.0[i][i] += C::new(c, 0.0);
    }
    out
}

pub fn cx() -> M4 {
    controlled(&x())
}
pub fn swap() -> M4 {
    let mut m = M4([[ZERO; 4]; 4]);
    for (i, j) in [(0, 0), (1, 2), (2, 1), (3, 3)] {
        m.0[i][j] = ONE;
    }
    m
}

/// How many qubits and parameters gate `name` takes, if it is a standard gate.
pub fn signature(name: &str) -> Option<(usize, usize)> {
    Some(match name {
        "id" | "x" | "y" | "z" | "h" | "s" | "sdg" | "t" | "tdg" | "sx" | "sxdg" => (1, 0),
        "rx" | "ry" | "rz" | "p" | "phase" | "u1" => (1, 1),
        "u2" => (1, 2),
        "u3" | "u" => (1, 3),
        "cx" | "cnot" | "cy" | "cz" | "ch" | "csx" | "swap" | "iswap" | "dcx" | "ecr" => (2, 0),
        "cp" | "cphase" | "cu1" | "crx" | "cry" | "crz" | "rxx" | "ryy" | "rzz" | "rzx" => (2, 1),
        "cu3" => (2, 3),
        "cu" => (2, 4),
        "ccx" | "cswap" => (3, 0),
        _ => return None,
    })
}

pub fn is_known_gate(name: &str) -> bool {
    signature(name).is_some()
}

/// The matrix of single-qubit gate `name`.
pub fn matrix_1q(name: &str, p: &[f64]) -> Option<M2> {
    Some(match (name, p) {
        ("id", []) => M2::IDENTITY,
        ("x", []) => x(),
        ("y", []) => y(),
        ("z", []) => z(),
        ("h", []) => h(),
        ("s", []) => s(),
        ("sdg", []) => s().adjoint(),
        ("t", []) => t(),
        ("tdg", []) => t().adjoint(),
        ("sx", []) => sx(),
        ("sxdg", []) => sx().adjoint(),
        ("rx", [a]) => rx(*a),
        ("ry", [a]) => ry(*a),
        ("rz", [a]) => rz(*a),
        ("p" | "phase" | "u1", [a]) => phase(*a),
        ("u2", [phi, lam]) => u3(FRAC_PI_2, *phi, *lam),
        ("u3" | "u", [a, b, c]) => u3(*a, *b, *c),
        _ => return None,
    })
}

/// The matrix of two-qubit gate `name`, big-endian (the first operand is the
/// most significant qubit, so a control comes first).
pub fn matrix_2q(name: &str, p: &[f64]) -> Option<M4> {
    Some(match (name, p) {
        ("cx" | "cnot", []) => cx(),
        ("cy", []) => controlled(&y()),
        ("cz", []) => controlled(&z()),
        ("ch", []) => controlled(&h()),
        ("csx", []) => controlled(&sx()),
        ("swap", []) => swap(),
        ("iswap", []) => {
            let mut m = M4::IDENTITY;
            m.0[1][1] = ZERO;
            m.0[2][2] = ZERO;
            m.0[1][2] = I;
            m.0[2][1] = I;
            m
        }
        ("dcx", []) => cx().swapped() * cx(),
        // Qiskit's definition: rzx(π/4) a, b; x a; rzx(-π/4) a, b.
        ("ecr", []) => {
            pauli_rotation(2, 0, -FRAC_PI_4)
                * kron(&x(), &M2::IDENTITY)
                * pauli_rotation(2, 0, FRAC_PI_4)
        }
        ("cp" | "cphase" | "cu1", [a]) => controlled(&phase(*a)),
        ("crx", [a]) => controlled(&rx(*a)),
        ("cry", [a]) => controlled(&ry(*a)),
        ("crz", [a]) => controlled(&rz(*a)),
        ("rxx", [a]) => pauli_rotation(0, 0, *a),
        ("ryy", [a]) => pauli_rotation(1, 1, *a),
        ("rzz", [a]) => pauli_rotation(2, 2, *a),
        ("rzx", [a]) => pauli_rotation(2, 0, *a),
        ("cu3", [a, b, c]) => controlled(&u3(*a, *b, *c)),
        ("cu", [a, b, c, g]) => controlled(&u3(*a, *b, *c).scale(C::from_polar(1.0, *g))),
        _ => return None,
    })
}
