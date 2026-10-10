"""Peephole optimizations."""

from __future__ import annotations

import math

import numpy as np

from qevedo.compiler.device import DeviceSpec, default_device
from qevedo.compiler.ir import Circuit, Instruction, Qubit
from qevedo.compiler.passes.base import Pass
from qevedo.compiler.synthesis.gates import gate_matrix, is_known_gate
from qevedo.compiler.synthesis.one_qubit import Basis1q, bases_for, synthesize_1q

__all__ = [
    "CancelAdjacentInverses",
    "CommutativeCancellation",
    "MergeSingleQubitGates",
    "merge_single_qubit_runs",
]

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


# How each gate acts on each of its qubits: "Z" if it commutes with Z there
# (diagonal on that qubit, like a CX control), "X" if it commutes with X (like
# a CX target). Two gates commute when every qubit they share has the same
# letter in both.
_Z_GATES = {
    "rz", "p", "u1", "phase", "z", "s", "sdg", "t", "tdg", "id",
    "cz", "cp", "cphase", "cu1", "crz", "rzz", "ccz",
}  # fmt: skip
_X_GATES = {"x", "sx", "sxdg", "rx", "rxx"}
_CONTROLLED_X = {"cx": 1, "cnot": 1, "crx": 1, "ccx": 2}
_CONTROLLED = {"cy": 1, "ch": 1, "cry": 1, "csx": 1, "cu": 1, "cu3": 1, "cswap": 1}
# Rotations that merge by adding their angles when they meet.
_MERGEABLE = {"rz", "rx", "ry", "p", "u1", "phase", "rzz", "rxx", "ryy", "cp", "cphase", "cu1"}


def _actions(inst: Instruction) -> list[str | None]:
    name = inst.name.lower()
    n = len(inst.qubits)
    if name in _Z_GATES:
        return ["Z"] * n
    if name in _X_GATES:
        return ["X"] * n
    if name in _CONTROLLED_X:
        controls = _CONTROLLED_X[name]
        return ["Z"] * controls + ["X"] * (n - controls)
    if name in _CONTROLLED:
        controls = _CONTROLLED[name]
        return ["Z"] * controls + [None] * (n - controls)
    return [None] * n


def _commute(a: Instruction, b: Instruction) -> bool:
    if a.clbits or b.clbits or a.name.lower() in NON_GATES or b.name.lower() in NON_GATES:
        return not set(a.qubits) & set(b.qubits)
    actions_a = dict(zip(a.qubits, _actions(a)))
    actions_b = dict(zip(b.qubits, _actions(b)))
    for q in set(a.qubits) & set(b.qubits):
        if actions_a[q] is None or actions_a[q] != actions_b[q]:
            return False
    return True


class CommutativeCancellation(Pass):
    """Cancel inverse pairs and merge rotations across gates they commute with.

    ``rz`` on a CX's control commutes with the CX, so ``rz(a) q0; cx q0, q1;
    rz(b) q0;`` becomes ``rz(a + b) q0; cx q0, q1;``, and two CXs with only
    gates that commute with them in between cancel.
    """

    name = "commutative_cancellation"

    #: How far back along a qubit to look for a partner.
    window = 32

    def run(self, circuit: Circuit, device: DeviceSpec | None = None) -> Circuit:
        del device
        out = circuit.copy()
        kept: list[Instruction | None] = []
        timelines: dict[Qubit, list[int]] = {}
        for inst in out.instructions:
            partner = self._partner(inst, kept, timelines)
            if partner is not None:
                previous = kept[partner]
                assert previous is not None
                merged = _combine(previous, inst)
                if merged is not _NO_MERGE:
                    if merged is None:
                        kept[partner] = None
                        for q in previous.qubits:
                            timelines[q].remove(partner)
                    else:
                        kept[partner] = merged
                    continue
            for q in inst.qubits:
                timelines.setdefault(q, []).append(len(kept))
            kept.append(inst)
        out.instructions = [inst for inst in kept if inst is not None]
        return out

    def _partner(
        self,
        inst: Instruction,
        kept: list[Instruction | None],
        timelines: dict[Qubit, list[int]],
    ) -> int | None:
        """An earlier gate on the same qubits that ``inst`` can reach by commuting."""
        if not inst.qubits or inst.clbits or inst.name.lower() in NON_GATES:
            return None
        first = timelines.get(inst.qubits[0], [])
        for index in reversed(first[-self.window :]):
            candidate = kept[index]
            assert candidate is not None
            if (
                set(candidate.qubits) == set(inst.qubits)
                and _combine(candidate, inst) is not _NO_MERGE
                and _reachable(inst, index, kept, timelines)
            ):
                return index
            if not _commute(candidate, inst):
                return None
        return None


def _reachable(
    inst: Instruction, index: int, kept: list[Instruction | None], timelines: dict[Qubit, list[int]]
) -> bool:
    """Whether every gate after ``index`` on ``inst``'s other qubits commutes with ``inst``."""
    return all(
        _commute(kept[j], inst)  # type: ignore[arg-type]
        for q in inst.qubits[1:]
        for j in timelines[q]
        if j > index
    )


_NO_MERGE = object()


def _combine(left: Instruction, right: Instruction):
    """``left`` followed by ``right`` as one gate, None if they cancel, or _NO_MERGE."""
    if _is_inverse_pair(left, right):
        return None
    a, b = left.name.lower(), right.name.lower()
    if a != b or a not in _MERGEABLE or len(left.params) != 1 or len(right.params) != 1:
        return _NO_MERGE
    if left.qubits != right.qubits and not (
        a in _SYMMETRIC and set(left.qubits) == set(right.qubits)
    ):
        return _NO_MERGE
    angle = left.params[0] + right.params[0]
    period = _ROTATIONS.get(a, 2 * math.pi)
    if abs(math.remainder(angle, period)) < 1e-12:
        return None
    return Instruction(left.name, left.qubits, (angle,))
