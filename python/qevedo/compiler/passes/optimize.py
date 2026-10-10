"""Peephole optimizations."""

from __future__ import annotations

import math

import numpy as np

from qevedo.compiler.device import DeviceSpec, default_device
from qevedo.compiler.ir import Circuit, Instruction, Qubit
from qevedo.compiler.passes.base import Pass
from qevedo.compiler.synthesis.gates import gate_matrix, is_known_gate
from qevedo.compiler.synthesis.one_qubit import Basis1q, bases_for, synthesize_1q

__all__ = ["CancelAdjacentInverses", "MergeSingleQubitGates", "merge_single_qubit_runs"]

NON_GATES = {"barrier", "measure", "reset"}

# Gates equal to their own inverse.
_SELF_INVERSE = {"x", "y", "z", "h", "id", "cx", "cy", "cz", "ch", "swap", "ecr", "ccx", "cswap"}
# Pairs of fixed gates that are each other's inverse.
_INVERSE = {"s": "sdg", "sdg": "s", "t": "tdg", "tdg": "t", "sx": "sxdg", "sxdg": "sx"}
# Rotations whose inverse is the same gate with the angle negated, with the
# period of the angle: a controlled rotation by 2π is a CZ-like phase, not a no-op.
_ROTATIONS = {
    **dict.fromkeys(
        ["rx", "ry", "rz", "p", "u1", "phase", "rxx", "ryy", "rzz", "rzx", "cp"], 2 * math.pi
    ),
    **dict.fromkeys(["crx", "cry", "crz"], 4 * math.pi),
}
# Two-qubit gates that do not depend on the order of their qubits.
_SYMMETRIC = {"cz", "swap", "rxx", "ryy", "rzz", "cp", "cphase"}


def _angles_close(a: float, b: float, period: float, tol: float = 1e-12) -> bool:
    return abs(math.remainder(a - b, period)) < tol


class CancelAdjacentInverses(Pass):
    """Remove pairs of gates that undo each other with nothing in between on their qubits."""

    name = "optimize"

    def run(self, circuit: Circuit, device: DeviceSpec | None = None) -> Circuit:
        del device
        out = circuit.copy()
        kept: list[Instruction | None] = []
        # For each qubit, the indices in `kept` of the instructions on it, in order.
        stacks: dict[Qubit, list[int]] = {}
        for inst in out.instructions:
            tops = {stacks[q][-1] if stacks.get(q) else None for q in inst.qubits}
            if len(tops) == 1 and (top := tops.pop()) is not None:
                previous = kept[top]
                assert previous is not None
                if set(previous.qubits) == set(inst.qubits) and _is_inverse_pair(previous, inst):
                    kept[top] = None
                    for q in inst.qubits:
                        stacks[q].pop()
                    continue
            for q in inst.qubits:
                stacks.setdefault(q, []).append(len(kept))
            kept.append(inst)
        out.instructions = [inst for inst in kept if inst is not None]
        return out


def _is_inverse_pair(left: Instruction, right: Instruction) -> bool:
    a, b = left.name.lower(), right.name.lower()
    if a in NON_GATES or b in NON_GATES or left.clbits or right.clbits:
        return False
    if left.qubits != right.qubits and not (a in _SYMMETRIC and a == b):
        return False
    if a == b and a in _SELF_INVERSE:
        return True
    if _INVERSE.get(a) == b:
        return True
    if a == b and a in _ROTATIONS and len(left.params) == 1 == len(right.params):
        return _angles_close(left.params[0], -right.params[0], _ROTATIONS[a])
    return False


class MergeSingleQubitGates(Pass):
    """Replace each run of single-qubit gates on a qubit by its shortest native equivalent."""

    name = "merge_1q"

    def run(self, circuit: Circuit, device: DeviceSpec | None = None) -> Circuit:
        native = {g.lower() for g in (device or default_device()).native_gates}
        out = circuit.copy()
        out.instructions = merge_single_qubit_runs(out.instructions, native)
        return out


def merge_single_qubit_runs(
    instructions: list[Instruction], native: set[str], bases: list[Basis1q] | None = None
) -> list[Instruction]:
    """``instructions`` with each run of single-qubit gates on a qubit resynthesized.

    A run is kept as it is when it is already native and no shorter. Runs on
    different qubits commute, so each is emitted where its qubit is next used.
    """
    if bases is None:
        bases = bases_for(native)
    result: list[Instruction] = []
    runs: dict[Qubit, list[Instruction]] = {}

    def flush(qubit: Qubit) -> None:
        run = runs.pop(qubit, None)
        if run:
            result.extend(_resynthesize(run, bases, native))

    for inst in instructions:
        name = inst.name.lower()
        if len(inst.qubits) == 1 and name not in NON_GATES and is_known_gate(name):
            runs.setdefault(inst.qubits[0], []).append(inst)
            continue
        for q in inst.qubits:
            flush(q)
        result.append(inst)
    for q in list(runs):
        flush(q)
    return result


def _resynthesize(
    run: list[Instruction], bases: list[Basis1q], native: set[str]
) -> list[Instruction]:
    if not bases:
        return run
    matrix = np.eye(2, dtype=complex)
    for inst in run:
        matrix = gate_matrix(inst.name.lower(), inst.params) @ matrix
    ops = synthesize_1q(matrix, bases)
    if all(inst.name.lower() in native for inst in run) and len(ops) >= len(run):
        return run
    qubit = run[0].qubits[0]
    return [Instruction(name, (qubit,), params) for name, params in ops]
