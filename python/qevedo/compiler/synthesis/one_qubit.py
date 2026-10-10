"""Optimal single-qubit synthesis into a device's native gates.

Every single-qubit unitary is, up to global phase, ``Rz(phi) Ry(theta) Rz(lam)``
(a ZYZ Euler decomposition). Each native basis below turns those three angles
into the fewest gates that basis allows, dropping rotations by multiples of
2π, so a gate that needs no pulse at all disappears.
"""

from __future__ import annotations

import cmath
import math
from collections.abc import Callable, Collection
from dataclasses import dataclass

import numpy as np

from qevedo.compiler.synthesis.gates import gate_matrix

__all__ = ["Basis1q", "Op", "bases_for", "euler_zyz", "synthesize_1q"]

#: One gate of a synthesized sequence: (name, params).
Op = tuple[str, tuple[float, ...]]

TOL = 1e-10


def _wrap(angle: float) -> float:
    """``angle`` in (-π, π]."""
    wrapped = math.remainder(angle, 2 * math.pi)
    return math.pi if math.isclose(wrapped, -math.pi, abs_tol=TOL) else wrapped


def _is_zero(angle: float) -> bool:
    return abs(_wrap(angle)) < TOL


def euler_zyz(u: np.ndarray) -> tuple[float, float, float]:
    """Angles ``(theta, phi, lam)`` with ``u ≅ Rz(phi) Ry(theta) Rz(lam)`` up to global phase.

    ``theta`` is in [0, π].
    """
    su = u / cmath.sqrt(np.linalg.det(u))
    # su = [[e^{-i(phi+lam)/2} cos(theta/2), -e^{-i(phi-lam)/2} sin(theta/2)],
    #       [e^{ i(phi-lam)/2} sin(theta/2),  e^{ i(phi+lam)/2} cos(theta/2)]]
    a, b = su[1, 1], su[1, 0]
    theta = 2 * math.atan2(abs(b), abs(a))
    plus = 2 * cmath.phase(a) if abs(a) > TOL else 0.0
    minus = 2 * cmath.phase(b) if abs(b) > TOL else 0.0
    return theta, (plus + minus) / 2, (plus - minus) / 2


@dataclass(frozen=True)
class Basis1q:
    """A family of native single-qubit gates and how to write a unitary in it."""

    name: str
    emit: Callable[[np.ndarray], list[Op]]


def _rz(z: str, angle: float) -> list[Op]:
    return [] if _is_zero(angle) else [(z, (_wrap(angle),))]


def _zsx(z: str, x_gate: str | None) -> Callable[[np.ndarray], list[Op]]:
    def emit(u: np.ndarray) -> list[Op]:
        theta, phi, lam = euler_zyz(u)
        if _is_zero(theta):
            return _rz(z, phi + lam)
        if abs(theta - math.pi / 2) < TOL:
            # U3(π/2, φ, λ) ≅ Rz(φ + π/2) SX Rz(λ - π/2)
            return _rz(z, lam - math.pi / 2) + [("sx", ())] + _rz(z, phi + math.pi / 2)
        if x_gate and abs(theta - math.pi) < TOL:
            # U3(π, φ, λ) ≅ Rz(φ + π) X Rz(λ) = Rz(φ - λ + π) X
            return [(x_gate, ())] + _rz(z, phi - lam + math.pi)
        # U3(θ, φ, λ) ≅ Rz(φ + π) SX Rz(θ + π) SX Rz(λ)
        return (
            _rz(z, lam)
            + [("sx", ())]
            + _rz(z, theta + math.pi)
            + [("sx", ())]
            + _rz(z, phi + math.pi)
        )

    return emit


def _zxz(z: str) -> Callable[[np.ndarray], list[Op]]:
    # Ry(θ) = Rz(π/2) Rx(θ) Rz(-π/2)
    def emit(u: np.ndarray) -> list[Op]:
        theta, phi, lam = euler_zyz(u)
        if _is_zero(theta):
            return _rz(z, phi + lam)
        return _rz(z, lam - math.pi / 2) + [("rx", (_wrap(theta),))] + _rz(z, phi + math.pi / 2)

    return emit


def _zyz(z: str) -> Callable[[np.ndarray], list[Op]]:
    def emit(u: np.ndarray) -> list[Op]:
        theta, phi, lam = euler_zyz(u)
        if _is_zero(theta):
            return _rz(z, phi + lam)
        return _rz(z, lam) + [("ry", (_wrap(theta),))] + _rz(z, phi)

    return emit


def _xyx(u: np.ndarray) -> list[Op]:
    # H Rz(a) H = Rx(a) and H Ry(b) H = Ry(-b), so the ZYZ angles of H U H give U as Rx Ry Rx.
    h = gate_matrix("h")
    theta, phi, lam = euler_zyz(h @ u @ h)
    if _is_zero(theta):
        return _rz("rx", phi + lam)
    return _rz("rx", lam) + [("ry", (_wrap(-theta),))] + _rz("rx", phi)


def _u3(name: str) -> Callable[[np.ndarray], list[Op]]:
    def emit(u: np.ndarray) -> list[Op]:
        theta, phi, lam = euler_zyz(u)
        if _is_zero(theta) and _is_zero(phi + lam):
            return []
        return [(name, (_wrap(theta), _wrap(phi), _wrap(lam)))]

    return emit


# Gates that act as Rz up to a global phase, in order of preference.
_Z_GATES = ("rz", "p", "u1", "phase")


def bases_for(native: Collection[str]) -> list[Basis1q]:
    """The single-qubit bases a native gate set supports."""
    native = {g.lower() for g in native}
    z = next((g for g in _Z_GATES if g in native), None)
    out: list[Basis1q] = []
    for name in ("u3", "u"):
        if name in native:
            out.append(Basis1q(name, _u3(name)))
            break
    if z and "sx" in native:
        out.append(Basis1q("zsx", _zsx(z, "x" if "x" in native else None)))
    if z and "rx" in native:
        out.append(Basis1q("zxz", _zxz(z)))
    if z and "ry" in native:
        out.append(Basis1q("zyz", _zyz(z)))
    if "rx" in native and "ry" in native:
        out.append(Basis1q("xyx", _xyx))
    return out


def synthesize_1q(u: np.ndarray, bases: list[Basis1q]) -> list[Op]:
    """The shortest gate sequence (in time order) over ``bases`` equal to ``u`` up to phase."""
    if not bases:
        raise ValueError("the native gate set has no universal single-qubit basis")
    return min((basis.emit(u) for basis in bases), key=len)


def ops_matrix(ops: list[Op]) -> np.ndarray:
    """The unitary of a time-ordered single-qubit sequence."""
    out = np.eye(2, dtype=complex)
    for name, params in ops:
        out = gate_matrix(name, params) @ out
    return out
