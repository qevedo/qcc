"""Lightweight peephole optimizations."""

from __future__ import annotations

import math
from typing import List, Optional

from qevedo.compiler.device import DeviceSpec
from qevedo.compiler.ir import Circuit, Instruction
from qevedo.compiler.passes.base import Pass


def _angles_close(a: float, b: float, tol: float = 1e-12) -> bool:
    return abs((a - b + math.pi) % (2 * math.pi) - math.pi) < tol


class CancelAdjacentInverses(Pass):
    name = "optimize"

    def run(self, circuit: Circuit, device: Optional[DeviceSpec] = None) -> Circuit:
        del device
        out = circuit.copy()
        out.instructions = self._cancel(out.instructions)
        return out

    def _cancel(self, instructions: List[Instruction]) -> List[Instruction]:
        stack: List[Instruction] = []
        for inst in instructions:
            if inst.name.lower() == "barrier":
                stack.append(inst.copy())
                continue
            if stack and self._is_inverse_pair(stack[-1], inst):
                stack.pop()
                continue
            stack.append(inst.copy())
        return stack

    def _is_inverse_pair(self, left: Instruction, right: Instruction) -> bool:
        if left.name.lower() != right.name.lower():
            return False
        if left.qubits != right.qubits:
            return False
        name = left.name.lower()
        if name == "x":
            return True
        if name == "sx":
            return True
        if name == "cx":
            return True
        if name == "rz" and left.params and right.params:
            return _angles_close(left.params[0], -right.params[0])
        return False
