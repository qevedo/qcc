"""Two-qubit block resynthesis.

A run of gates that only touches one pair of qubits is a single 4×4 unitary,
however many gates it has. This pass collects maximal such runs ("blocks"),
multiplies each into its unitary and replaces it by a fresh KAK synthesis when
that has fewer two-qubit gates, or as many and fewer gates in total. A block of
two CXs with rotations between them that happen to make a local gate, for
example, collapses to single-qubit gates.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from qevedo.compiler.device import DeviceSpec, default_device
from qevedo.compiler.ir import Circuit, Instruction, Qubit
from qevedo.compiler.passes.base import Pass
from qevedo.compiler.passes.decompose import Lowering, cost
from qevedo.compiler.synthesis.gates import embed, gate_matrix, is_known_gate

__all__ = ["ResynthesizeTwoQubitBlocks"]

NON_GATES = {"barrier", "measure", "reset"}


@dataclass
class _Block:
    pair: tuple[Qubit, Qubit]
    instructions: list[Instruction] = field(default_factory=list)

    def matrix(self) -> np.ndarray:
        """The block's unitary, with ``pair[0]`` as the most significant qubit."""
        out = np.eye(4, dtype=complex)
        for inst in self.instructions:
            positions = tuple(self.pair.index(q) for q in inst.qubits)
            out = embed(gate_matrix(inst.name.lower(), inst.params), positions, 2) @ out
        return out


class ResynthesizeTwoQubitBlocks(Pass):
    """Replace each maximal two-qubit block by its cheapest synthesis."""

    name = "resynthesize_2q"

    def __init__(self) -> None:
        # Shared by the rounds of a pipeline, so blocks seen before are free.
        self._lowerings: dict[frozenset[str], Lowering] = {}

    def run(self, circuit: Circuit, device: DeviceSpec | None = None) -> Circuit:
        native = frozenset(g.lower() for g in (device or default_device()).native_gates)
        lowering = self._lowerings.setdefault(native, Lowering(native))
        out = circuit.copy()
        # The output, with a block standing where its first gate was: gates that
        # come later on other qubits commute with the rest of the block.
        result: list[Instruction | _Block] = []
        open_blocks: dict[Qubit, _Block] = {}
        # Single-qubit gates on a qubit that is in no block yet; a block that
        # starts on the qubit absorbs them.
        pending: dict[Qubit, list[Instruction]] = {}

        def close(qubit: Qubit) -> None:
            block = open_blocks.pop(qubit, None)
            if block is not None:
                for q in block.pair:
                    open_blocks.pop(q, None)
            result.extend(pending.pop(qubit, []))

        for inst in out.instructions:
            name = inst.name.lower()
            gate = name not in NON_GATES and not inst.clbits and is_known_gate(name)
            if gate and len(inst.qubits) == 1:
                (q,) = inst.qubits
                if q in open_blocks:
                    open_blocks[q].instructions.append(inst)
                else:
                    pending.setdefault(q, []).append(inst)
                continue
            if gate and len(inst.qubits) == 2:
                a, b = inst.qubits
                block = open_blocks.get(a)
                if block is not None and block is open_blocks.get(b):
                    block.instructions.append(inst)
                    continue
                close(a)
                close(b)
                block = _Block((a, b))
                block.instructions = pending.pop(a, []) + pending.pop(b, []) + [inst]
                open_blocks[a] = open_blocks[b] = block
                result.append(block)
                continue
            for q in inst.qubits:
                close(q)
            result.append(inst)
        for q in list(pending):
            result.extend(pending.pop(q))

        out.instructions = []
        for item in result:
            if isinstance(item, _Block):
                out.instructions.extend(self._best(item, lowering))
            else:
                out.instructions.append(item)
        return out

    @staticmethod
    def _best(block: _Block, lowering: Lowering) -> list[Instruction]:
        two_qubit_gates = sum(1 for inst in block.instructions if len(inst.qubits) == 2)
        if two_qubit_gates < 2 and len(block.instructions) <= 3:
            return block.instructions
        candidate = lowering.two_qubit(block.matrix(), block.pair)
        return candidate if cost(candidate) < cost(block.instructions) else block.instructions
