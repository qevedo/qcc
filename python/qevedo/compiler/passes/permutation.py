"""SWAP elision: remove SWAP gates by relabelling the qubits of what follows.

The reference for the Rust core's ``core/src/permutation.rs``, and what
Qiskit's ``ElidePermutations`` does. A SWAP costs three CNOTs; relabelling
costs nothing, but leaves the qubits' states on other wires at the end. That
cannot be seen when every moved state is measured last, so the pass then
always elides; otherwise it does only when asked to, and records where each
state ends in the circuit's ``output_permutation``.

Relabelled gates may act on new pairs of qubits, so the pass leaves circuits
for devices that are not all-to-all alone.
"""

from __future__ import annotations

from qevedo.compiler.device import DeviceSpec
from qevedo.compiler.ir import Circuit, Instruction, Qubit
from qevedo.compiler.passes.base import Pass

__all__ = ["ElideSwaps"]


class ElideSwaps(Pass):
    """Remove SWAP gates by relabelling.

    With ``relabel_outputs`` false (the default), only when the relabelling
    cannot be seen: every qubit whose state moves is measured last.
    """

    name = "elide_swaps"

    def __init__(self, relabel_outputs: bool = False):
        self.relabel_outputs = relabel_outputs

    def run(self, circuit: Circuit, device: DeviceSpec | None = None) -> Circuit:
        if device is not None and device.topology != "all_to_all":
            return circuit
        if not any(inst.name.lower() == "swap" for inst in circuit.instructions):
            return circuit
        # The wire that holds each qubit's state, for the qubits that moved.
        wires: dict[Qubit, Qubit] = dict(circuit.output_permutation)
        instructions = []
        for inst in circuit.instructions:
            qs = tuple(wires.get(q, q) for q in inst.qubits)
            if inst.name.lower() == "swap" and len(qs) == 2 and not inst.clbits:
                a, b = inst.qubits
                wires[a], wires[b] = qs[1], qs[0]
            else:
                instructions.append(
                    Instruction(inst.name, qs, inst.params, inst.clbits)
                )
        moved = {q: w for q, w in wires.items() if q != w}
        if not self.relabel_outputs and not _measured_last(
            instructions, set(moved.values())
        ):
            return circuit
        out = circuit.copy()
        out.instructions = instructions
        out.output_permutation = moved
        return out


def _measured_last(instructions: list[Instruction], wires: set[Qubit]) -> bool:
    """Whether every wire of ``wires`` is measured last (barriers aside)."""
    last = {}
    for inst in instructions:
        if inst.name.lower() != "barrier":
            for q in inst.qubits:
                last[q] = inst.name.lower()
    return all(last.get(w) == "measure" for w in wires)
