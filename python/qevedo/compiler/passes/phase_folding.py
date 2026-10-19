"""Phase folding: merge rotations that act on the same parity, however far apart.

From Amy, Maslov and Mosca (2014) and Nam, Ross, Su, Childs and Maslov (2018),
in the path-sum form of Amy's Feynman; the reference for the Rust core's
``core/src/phase_folding.rs``.

Every wire carries an affine parity: the XOR of some variables, plus a
constant. The qubits start as variables; a CNOT XORs its control's parity into
its target's, an X flips the constant, and any other gate that is not diagonal
gives its qubits fresh variables. A diagonal rotation (Rz, phase, T, S, Z, or
Rzz on two wires) multiplies the state by ``e^{iθ f}`` for the parity ``f`` it
acts on, and the circuit's phase is the sum of all such terms, so rotations on
the same parity merge into one: the last adjustable one takes the sum and the
others go. Rotations on a constant parity are global phases and go too.
"""

from __future__ import annotations

import math

from qevedo.compiler.device import DeviceSpec
from qevedo.compiler.ir import Circuit, Instruction
from qevedo.compiler.passes.base import Pass

__all__ = ["PhaseFolding"]

# The phase of diagonal single-qubit gates on their wire's parity, and whether
# their angle is a parameter that can be changed.
_FIXED = {"z": math.pi, "s": math.pi / 2, "sdg": -math.pi / 2, "t": math.pi / 4, "tdg": -math.pi / 4}
_ADJUSTABLE = {"rz", "p", "phase", "u1"}
# Diagonal gates that are not a phase on one parity: they leave parities alone.
_OTHER_DIAGONAL = {"cz", "cp", "cphase", "cu1", "crz", "id", "barrier"}


class PhaseFolding(Pass):
    """Merge rotations of the circuit that act on the same parity."""

    name = "phase_folding"

    def run(self, circuit: Circuit, device: DeviceSpec | None = None) -> Circuit:
        del device
        index = {}
        for inst in circuit.instructions:
            for q in inst.qubits:
                index.setdefault(q, len(index))
        # A parity is (bitset of variables as an int, constant bit).
        wires = {q: (1 << i, 0) for q, i in index.items()}
        next_variable = len(index)
        # Terms by parity: (instruction position, sign, coefficient, adjustable).
        groups: dict[int, list[tuple[int, float, float, bool]]] = {}

        def add(parity: tuple[int, int], position: int, theta: float, adjustable: bool) -> None:
            sign = -1.0 if parity[1] else 1.0
            groups.setdefault(parity[0], []).append((position, sign, sign * theta, adjustable))

        for position, inst in enumerate(circuit.instructions):
            name, qs = inst.name.lower(), inst.qubits
            gate = not inst.clbits and name not in ("measure", "reset")
            if gate and name in ("cx", "cnot") and len(qs) == 2:
                (vc, cc), (vt, ct) = wires[qs[0]], wires[qs[1]]
                wires[qs[1]] = (vt ^ vc, ct ^ cc)
            elif gate and name == "x":
                v, c = wires[qs[0]]
                wires[qs[0]] = (v, c ^ 1)
            elif gate and name == "swap":
                wires[qs[0]], wires[qs[1]] = wires[qs[1]], wires[qs[0]]
            elif gate and name == "rzz":
                (va, ca), (vb, cb) = wires[qs[0]], wires[qs[1]]
                add((va ^ vb, ca ^ cb), position, float(inst.params[0]), True)
            elif gate and name in _OTHER_DIAGONAL:
                pass
            elif gate and name in _ADJUSTABLE and len(qs) == 1:
                add(wires[qs[0]], position, float(inst.params[0]), True)
            elif gate and name in _FIXED:
                add(wires[qs[0]], position, _FIXED[name], False)
            else:
                for q in qs:
                    wires[q] = (1 << next_variable, 0)
                    next_variable += 1

        dropped: set[int] = set()
        angles: dict[int, float] = {}
        for variables, terms in groups.items():
            if variables == 0:  # a global phase
                dropped.update(t[0] for t in terms)
                continue
            keepers = [t for t in terms if t[3]]
            if len(terms) < 2 or not keepers:
                continue
            keeper = keepers[-1]
            dropped.update(t[0] for t in terms if t[0] != keeper[0])
            total = sum(t[2] for t in terms)
            angle = math.remainder(keeper[1] * total, 2 * math.pi)
            if abs(angle) < 1e-12:
                dropped.add(keeper[0])
            else:
                angles[keeper[0]] = angle

        out = circuit.copy()
        out.instructions = []
        for position, inst in enumerate(circuit.instructions):
            if position in dropped:
                continue
            if position in angles:
                inst = Instruction(inst.name, inst.qubits, (angles[position],), inst.clbits)
            out.instructions.append(inst)
        return out
