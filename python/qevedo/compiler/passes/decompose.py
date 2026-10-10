"""Lower every gate to the device's native gates, with the fewest gates.

Each gate is lowered in every way that applies, and the cheapest result (fewest
two-qubit gates, then fewest gates) is kept:

* single-qubit gates: an Euler decomposition into the device's single-qubit
  basis, which is optimal;
* two-qubit gates: a KAK decomposition into the device's two-qubit gate, with
  the optimal number of two-qubit gates (see :mod:`qevedo.compiler.synthesis`);
* gates of the standard libraries: their textbook definition, lowered
  recursively. These often need fewer single-qubit gates than a KAK
  decomposition (``crz`` is ``rz, cx, rz, cx``), and they are the only route
  for gates on three qubits.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Sequence
from functools import lru_cache

import numpy as np

from qevedo.compiler.device import DeviceSpec, default_device
from qevedo.compiler.ir import Circuit, Instruction, Qubit
from qevedo.compiler.passes.base import Pass
from qevedo.compiler.passes.optimize import merge_single_qubit_runs
from qevedo.compiler.synthesis.gates import gate_matrix, is_known_gate
from qevedo.compiler.synthesis.one_qubit import bases_for, count_1q, synthesize_1q
from qevedo.compiler.synthesis.two_qubit import best_2q_bases, synthesize_2q

__all__ = ["DecomposeToNative", "Lowering", "LoweringError"]

PI = math.pi

#: Operations that are not gates and pass through unchanged.
NON_GATES = {"barrier", "measure", "reset"}

Body = list[tuple[str, tuple[float, ...], tuple[int, ...]]]


def _crz(t: float) -> Body:
    return [("rz", (t / 2,), (1,)), ("cx", (), (0, 1)), ("rz", (-t / 2,), (1,)), ("cx", (), (0, 1))]


def _cp(t: float) -> Body:
    return [
        ("p", (t / 2,), (0,)),
        ("cx", (), (0, 1)),
        ("p", (-t / 2,), (1,)),
        ("cx", (), (0, 1)),
        ("p", (t / 2,), (1,)),
    ]


def _cu(theta: float, phi: float, lam: float, gamma: float = 0.0) -> Body:
    return [
        ("p", (gamma + (lam + phi) / 2,), (0,)),
        ("p", ((lam - phi) / 2,), (1,)),
        ("cx", (), (0, 1)),
        ("u3", (-theta / 2, 0.0, -(phi + lam) / 2), (1,)),
        ("cx", (), (0, 1)),
        ("u3", (theta / 2, phi, 0.0), (1,)),
    ]


def _rzz(t: float) -> Body:
    return [("cx", (), (0, 1)), ("rz", (t,), (1,)), ("cx", (), (0, 1))]


# Definitions of the standard gates (qelib1.inc, stdgates.inc and Qiskit's
# library) in smaller gates, as (name, params, operand indices) in time order.
DEFINITIONS: dict[str, Callable[..., Body]] = {
    "swap": lambda: [("cx", (), (0, 1)), ("cx", (), (1, 0)), ("cx", (), (0, 1))],
    "cz": lambda: [("h", (), (1,)), ("cx", (), (0, 1)), ("h", (), (1,))],
    "cy": lambda: [("sdg", (), (1,)), ("cx", (), (0, 1)), ("s", (), (1,))],
    "ch": lambda: [("ry", (PI / 4,), (1,)), ("cx", (), (0, 1)), ("ry", (-PI / 4,), (1,))],
    "csx": lambda: [("h", (), (1,))] + _cp(PI / 2) + [("h", (), (1,))],
    "crz": _crz,
    "crx": lambda t: [("h", (), (1,))] + _crz(t) + [("h", (), (1,))],
    "cry": lambda t: [
        ("ry", (t / 2,), (1,)),
        ("cx", (), (0, 1)),
        ("ry", (-t / 2,), (1,)),
        ("cx", (), (0, 1)),
    ],
    "cp": _cp,
    "cphase": _cp,
    "cu1": _cp,
    "cu3": _cu,
    "cu": _cu,
    "rzz": _rzz,
    "rxx": lambda t: (
        [("h", (), (0,)), ("h", (), (1,))] + _rzz(t) + [("h", (), (0,)), ("h", (), (1,))]
    ),
    "ryy": lambda t: (
        [("rx", (PI / 2,), (0,)), ("rx", (PI / 2,), (1,))]
        + _rzz(t)
        + [("rx", (-PI / 2,), (0,)), ("rx", (-PI / 2,), (1,))]
    ),
    "rzx": lambda t: [("h", (), (1,))] + _rzz(t) + [("h", (), (1,))],
    "dcx": lambda: [("cx", (), (0, 1)), ("cx", (), (1, 0))],
    "iswap": lambda: [
        ("s", (), (0,)),
        ("s", (), (1,)),
        ("h", (), (0,)),
        ("cx", (), (0, 1)),
        ("cx", (), (1, 0)),
        ("h", (), (1,)),
    ],
    # Toffoli with 6 CX, which is optimal.
    "ccx": lambda: [
        ("h", (), (2,)),
        ("cx", (), (1, 2)),
        ("tdg", (), (2,)),
        ("cx", (), (0, 2)),
        ("t", (), (2,)),
        ("cx", (), (1, 2)),
        ("tdg", (), (2,)),
        ("cx", (), (0, 2)),
        ("t", (), (1,)),
        ("t", (), (2,)),
        ("h", (), (2,)),
        ("cx", (), (0, 1)),
        ("t", (), (0,)),
        ("tdg", (), (1,)),
        ("cx", (), (0, 1)),
    ],
    "cswap": lambda: [("cx", (), (2, 1)), ("ccx", (), (0, 1, 2)), ("cx", (), (2, 1))],
    # Toffoli up to a relative phase, with 3 CX.
    "rccx": lambda: [
        ("u2", (0.0, PI), (2,)),
        ("u1", (PI / 4,), (2,)),
        ("cx", (), (1, 2)),
        ("u1", (-PI / 4,), (2,)),
        ("cx", (), (0, 2)),
        ("u1", (PI / 4,), (2,)),
        ("cx", (), (1, 2)),
        ("u1", (-PI / 4,), (2,)),
        ("u2", (0.0, PI), (2,)),
    ],
}


class LoweringError(ValueError):
    """A gate cannot be written in the device's native gates."""


class Lowering:
    """Lowers single instructions to a fixed native gate set."""

    def __init__(self, native_gates: Sequence[str] | set[str]):
        self.native = {g.lower() for g in native_gates}
        self.bases_1q = bases_for(self.native)
        self.bases_2q = best_2q_bases(self.native)
        self._cost_1q = lru_cache(maxsize=None)(self._cost_1q_uncached)
        self._memo: dict[tuple[str, tuple[float, ...], int], list[Instruction]] = {}
        self._memo_2q: dict[bytes, list[Instruction]] = {}

    def lower(self, inst: Instruction) -> list[Instruction]:
        name = inst.name.lower()
        if name in NON_GATES or name in self.native:
            return [inst.copy()]
        # Lower each distinct gate once, on placeholder qubits, then relabel.
        key = (name, tuple(inst.params), len(inst.qubits))
        if key not in self._memo:
            placeholders = tuple(Qubit("", i) for i in range(len(inst.qubits)))
            self._memo[key] = self._lower(Instruction(name, placeholders, inst.params))
        return [
            Instruction(low.name, tuple(inst.qubits[q.index] for q in low.qubits), low.params)
            for low in self._memo[key]
        ]

    def _lower(self, inst: Instruction) -> list[Instruction]:
        name = inst.name
        candidates: list[list[Instruction]] = []
        if name in DEFINITIONS:
            out: list[Instruction] = []
            for sub, params, operands in DEFINITIONS[name](*inst.params):
                out.extend(
                    self.lower(Instruction(sub, tuple(inst.qubits[i] for i in operands), params))
                )
            candidates.append(merge_single_qubit_runs(out, self.native, self.bases_1q))
        if is_known_gate(name) and len(inst.qubits) <= 2:
            matrix = gate_matrix(name, inst.params)
            if len(inst.qubits) == 1:
                candidates.append(self._one_qubit(matrix, inst.qubits[0]))
            elif self.bases_2q:
                candidates.append(self.two_qubit(matrix, inst.qubits))
        if not candidates:
            if not is_known_gate(name) and name not in DEFINITIONS:
                raise LoweringError(
                    f"cannot lower gate {inst.name!r}: qcc does not know its definition yet"
                )
            raise LoweringError(
                f"cannot lower {inst.name!r}: the native gates {sorted(self.native)} include no "
                "supported two-qubit gate (cx, cz, cy, ch, ecr, rxx, ryy, rzz or rzx)"
            )
        return min(candidates, key=cost)

    def _one_qubit(self, matrix: np.ndarray, qubit: Qubit) -> list[Instruction]:
        if not self.bases_1q:
            raise LoweringError(
                f"the native gates {sorted(self.native)} have no universal single-qubit basis"
            )
        ops = synthesize_1q(matrix, self.bases_1q)
        return [Instruction(name, (qubit,), params) for name, params in ops]

    def _cost_1q_uncached(self, key: bytes) -> int:
        matrix = np.frombuffer(key, dtype=complex).reshape(2, 2)
        return count_1q(matrix, self.bases_1q) if self.bases_1q else 0

    def two_qubit(self, matrix: np.ndarray, qubits: Sequence[Qubit]) -> list[Instruction]:
        """The cheapest native sequence for the two-qubit unitary ``matrix`` on ``qubits``."""
        key = np.ascontiguousarray(matrix, dtype=complex).tobytes()
        if key not in self._memo_2q:
            placeholders = (Qubit("", 0), Qubit("", 1))
            self._memo_2q[key] = self._two_qubit(matrix, placeholders)
        return [
            Instruction(low.name, tuple(qubits[q.index] for q in low.qubits), low.params)
            for low in self._memo_2q[key]
        ]

    def _two_qubit(self, matrix: np.ndarray, qubits: Sequence[Qubit]) -> list[Instruction]:

        def cost_1q(u: np.ndarray) -> int:
            return self._cost_1q(np.ascontiguousarray(u, dtype=complex).tobytes())

        best: list[Instruction] | None = None
        for basis in self.bases_2q:
            out: list[Instruction] = []
            for name, params, operands in synthesize_2q(matrix, basis, cost_1q):
                if name == "u":
                    out.extend(self._one_qubit(params, qubits[operands[0]]))
                else:
                    out.append(Instruction(name, tuple(qubits[i] for i in operands), params))
            if best is None or cost(out) < cost(best):
                best = out
        assert best is not None
        return best


def cost(instructions: Sequence[Instruction]) -> tuple[int, int]:
    """Two-qubit gates first, then all gates."""
    return (sum(1 for inst in instructions if len(inst.qubits) > 1), len(instructions))


class DecomposeToNative(Pass):
    """Rewrite every gate that is not native to the device in native gates."""

    name = "decompose"

    def run(self, circuit: Circuit, device: DeviceSpec | None = None) -> Circuit:
        lowering = Lowering((device or default_device()).native_gates)
        out = circuit.copy()
        out.instructions = []
        for inst in circuit.instructions:
            out.instructions.extend(lowering.lower(inst))
        return out
