"""Decompose abstract gates into a native {rz, sx, x, cx} basis."""

from __future__ import annotations

import math
from typing import List, Optional

from qcc.device import DeviceSpec
from qcc.ir import Circuit, Instruction, Qubit
from qcc.passes.base import Pass

PI = math.pi
HALF_PI = PI / 2


def _rz(q: Qubit, angle: float) -> Instruction:
    return Instruction("rz", (q,), (angle,))


def _sx(q: Qubit) -> Instruction:
    return Instruction("sx", (q,))


def _x(q: Qubit) -> Instruction:
    return Instruction("x", (q,))


def decompose_u(q: Qubit, theta: float, phi: float, lam: float) -> List[Instruction]:
    """Decompose U(θ, φ, λ) into rz/sx gates (global phase ignored)."""
    return [
        _rz(q, lam),
        _sx(q),
        _rz(q, theta),
        _sx(q),
        _rz(q, phi),
    ]


def decompose_h(q: Qubit) -> List[Instruction]:
    """Decompose Hadamard into rz/sx."""
    return [_rz(q, -HALF_PI), _sx(q), _rz(q, -HALF_PI)]


class DecomposeToNative(Pass):
    name = "decompose"

    def run(self, circuit: Circuit, device: Optional[DeviceSpec] = None) -> Circuit:
        del device
        out = circuit.copy()
        out.instructions = []
        for inst in circuit.instructions:
            out.instructions.extend(self._decompose(inst))
        return out

    def _decompose(self, inst: Instruction) -> List[Instruction]:
        name = inst.name.lower()
        if name in {"barrier", "measure", "reset", "cx", "rz", "sx", "x", "id"}:
            return [inst.copy()]
        if name in {"h", "ch"}:
            return decompose_h(inst.qubits[0])
        if name in {"u", "u3"}:
            if len(inst.params) != 3:
                raise ValueError(f"{inst.name} expects 3 parameters, got {inst.params}")
            return decompose_u(inst.qubits[0], *inst.params)
        if name == "u2":
            theta, phi = (0.0, 0.0) if not inst.params else (inst.params[0], inst.params[1])
            return decompose_u(inst.qubits[0], HALF_PI, phi, theta)
        if name == "u1":
            lam = inst.params[0] if inst.params else 0.0
            return [_rz(inst.qubits[0], lam)]
        if name == "p":
            return [_rz(inst.qubits[0], inst.params[0] if inst.params else 0.0)]
        if name == "z":
            return [_rz(inst.qubits[0], PI)]
        if name == "s":
            return [_rz(inst.qubits[0], HALF_PI)]
        if name == "t":
            return [_rz(inst.qubits[0], PI / 4)]
        if name == "y":
            return [_rz(inst.qubits[0], PI), _x(inst.qubits[0]), _rz(inst.qubits[0], PI)]
        # Pass through unknown gates; later passes or the backend may reject them.
        return [inst.copy()]
