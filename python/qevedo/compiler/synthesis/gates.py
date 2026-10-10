"""Unitary matrices of the standard OpenQASM gates.

Matrices are big-endian: for a gate applied to qubits ``(a, b, ...)``, ``a``
is the most significant bit of the row and column index, so a controlled gate
has its control first. The definitions follow ``qelib1.inc`` and
``stdgates.inc`` (and Qiskit for ``ecr``, ``iswap``, ``rzx``), including their
global phases.
"""

from __future__ import annotations

import cmath
import math
from collections.abc import Callable, Sequence

import numpy as np

__all__ = ["GATE_ARITY", "gate_matrix", "is_known_gate", "rx", "ry", "rz", "u3"]

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
H = np.array([[1, 1], [1, -1]], dtype=complex) / math.sqrt(2)
S = np.diag([1, 1j])
T = np.diag([1, cmath.exp(1j * math.pi / 4)])
SX = np.array([[1 + 1j, 1 - 1j], [1 - 1j, 1 + 1j]]) / 2


def rx(theta: float) -> np.ndarray:
    c, s = math.cos(theta / 2), math.sin(theta / 2)
    return np.array([[c, -1j * s], [-1j * s, c]])


def ry(theta: float) -> np.ndarray:
    c, s = math.cos(theta / 2), math.sin(theta / 2)
    return np.array([[c, -s], [s, c]], dtype=complex)


def rz(theta: float) -> np.ndarray:
    return np.diag([cmath.exp(-0.5j * theta), cmath.exp(0.5j * theta)])


def phase(lam: float) -> np.ndarray:
    return np.diag([1, cmath.exp(1j * lam)])


def u3(theta: float, phi: float, lam: float) -> np.ndarray:
    c, s = math.cos(theta / 2), math.sin(theta / 2)
    return np.array(
        [
            [c, -cmath.exp(1j * lam) * s],
            [cmath.exp(1j * phi) * s, cmath.exp(1j * (phi + lam)) * c],
        ]
    )


def controlled(u: np.ndarray) -> np.ndarray:
    """``u`` controlled on one extra qubit, the control being first."""
    n = u.shape[0]
    out = np.eye(2 * n, dtype=complex)
    out[n:, n:] = u
    return out


def _pauli_rotation(p: np.ndarray, q: np.ndarray, theta: float) -> np.ndarray:
    """exp(-i theta/2 P⊗Q)."""
    pq = np.kron(p, q)
    return math.cos(theta / 2) * np.eye(4) - 1j * math.sin(theta / 2) * pq


def _circuit(n: int, ops: Sequence[tuple[np.ndarray, Sequence[int]]]) -> np.ndarray:
    """Unitary of a small circuit given as (matrix, qubits) in time order."""
    out = np.eye(2**n, dtype=complex)
    for matrix, qubits in ops:
        out = embed(matrix, qubits, n) @ out
    return out


def embed(matrix: np.ndarray, qubits: Sequence[int], n: int) -> np.ndarray:
    """``matrix`` acting on ``qubits`` of an ``n``-qubit register (big-endian)."""
    k = len(qubits)
    rest = [q for q in range(n) if q not in qubits]
    order = list(qubits) + rest
    full = np.kron(matrix, np.eye(2 ** (n - k))).reshape([2] * (2 * n))
    # Axes of `full` are (order..., order...); move them back to 0..n-1.
    inverse = [order.index(q) for q in range(n)]
    full = full.transpose(inverse + [n + i for i in inverse])
    return full.reshape(2**n, 2**n)


SWAP = np.eye(4, dtype=complex)[[0, 2, 1, 3]]
CX = controlled(X)

_RZX_PLUS = _pauli_rotation(Z, X, math.pi / 4)
_RZX_MINUS = _pauli_rotation(Z, X, -math.pi / 4)

_FIXED: dict[str, np.ndarray] = {
    "id": I2,
    "x": X,
    "y": Y,
    "z": Z,
    "h": H,
    "s": S,
    "sdg": S.conj().T,
    "t": T,
    "tdg": T.conj().T,
    "sx": SX,
    "sxdg": SX.conj().T,
    "cx": CX,
    "cnot": CX,
    "cy": controlled(Y),
    "cz": controlled(Z),
    "ch": controlled(H),
    "csx": controlled(SX),
    "swap": SWAP,
    "iswap": np.array([[1, 0, 0, 0], [0, 0, 1j, 0], [0, 1j, 0, 0], [0, 0, 0, 1]]),
    "dcx": _circuit(2, [(CX, (0, 1)), (CX, (1, 0))]),
    # Qiskit's definition: rzx(pi/4) a, b; x a; rzx(-pi/4) a, b;
    "ecr": _RZX_MINUS @ np.kron(X, I2) @ _RZX_PLUS,
    "ccx": controlled(controlled(X)),
    "cswap": controlled(SWAP),
}

_PARAMETRIC: dict[str, tuple[int, Callable[..., np.ndarray]]] = {
    "rx": (1, rx),
    "ry": (1, ry),
    "rz": (1, rz),
    "p": (1, phase),
    "phase": (1, phase),
    "u1": (1, phase),
    "u2": (2, lambda phi, lam: u3(math.pi / 2, phi, lam)),
    "u3": (3, u3),
    "u": (3, u3),
    "cp": (1, lambda lam: controlled(phase(lam))),
    "cphase": (1, lambda lam: controlled(phase(lam))),
    "cu1": (1, lambda lam: controlled(phase(lam))),
    "crx": (1, lambda t: controlled(rx(t))),
    "cry": (1, lambda t: controlled(ry(t))),
    "crz": (1, lambda t: controlled(rz(t))),
    "cu3": (3, lambda t, p, lam: controlled(u3(t, p, lam))),
    "cu": (4, lambda t, p, lam, g: controlled(cmath.exp(1j * g) * u3(t, p, lam))),
    "rxx": (1, lambda t: _pauli_rotation(X, X, t)),
    "ryy": (1, lambda t: _pauli_rotation(Y, Y, t)),
    "rzz": (1, lambda t: _pauli_rotation(Z, Z, t)),
    "rzx": (1, lambda t: _pauli_rotation(Z, X, t)),
}

GATE_ARITY: dict[str, int] = {name: int(math.log2(m.shape[0])) for name, m in _FIXED.items()}
GATE_ARITY.update(
    {name: int(math.log2(f(*[0.0] * k).shape[0])) for name, (k, f) in _PARAMETRIC.items()}
)


def is_known_gate(name: str) -> bool:
    return name in _FIXED or name in _PARAMETRIC


def gate_matrix(name: str, params: Sequence[float] = ()) -> np.ndarray:
    """The unitary of gate ``name`` with ``params``; raises ``KeyError`` if unknown."""
    if name in _FIXED:
        if params:
            raise ValueError(f"{name} takes no parameters, got {len(params)}")
        return _FIXED[name]
    count, make = _PARAMETRIC[name]
    if len(params) != count:
        raise ValueError(f"{name} takes {count} parameter(s), got {len(params)}")
    return make(*params)
