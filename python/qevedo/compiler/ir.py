"""QCC intermediate representation."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple


@dataclass(frozen=True)
class Qubit:
    """A logical qubit identified by register name and index."""

    register: str
    index: int

    def label(self) -> str:
        return f"{self.register}[{self.index}]"


@dataclass(frozen=True)
class Clbit:
    """A classical bit identified by register name and index."""

    register: str
    index: int

    def label(self) -> str:
        return f"{self.register}[{self.index}]"


@dataclass
class Instruction:
    """One operation in a quantum circuit."""

    name: str
    qubits: Tuple[Qubit, ...] = ()
    params: Tuple[float, ...] = ()
    clbits: Tuple[Clbit, ...] = ()
    # Optional physical mapping filled in by layout passes.
    layout: Optional[Tuple[int, ...]] = None

    def copy(self) -> Instruction:
        return Instruction(
            name=self.name,
            qubits=self.qubits,
            params=self.params,
            clbits=self.clbits,
            layout=self.layout,
        )

    def is_two_qubit(self) -> bool:
        return len(self.qubits) == 2 and self.name.lower() not in {"measure", "barrier"}


@dataclass
class Circuit:
    """Ordered quantum + classical program."""

    num_qubits: int = 0
    num_clbits: int = 0
    qreg_sizes: Dict[str, int] = field(default_factory=dict)
    creg_sizes: Dict[str, int] = field(default_factory=dict)
    instructions: List[Instruction] = field(default_factory=list)

    def append(self, inst: Instruction) -> None:
        self.instructions.append(inst)

    def extend(self, insts: Sequence[Instruction]) -> None:
        self.instructions.extend(insts)

    def copy(self) -> Circuit:
        return Circuit(
            num_qubits=self.num_qubits,
            num_clbits=self.num_clbits,
            qreg_sizes=dict(self.qreg_sizes),
            creg_sizes=dict(self.creg_sizes),
            instructions=[inst.copy() for inst in self.instructions],
        )

    def gate_counts(self) -> Dict[str, int]:
        counts: Dict[str, int] = {}
        for inst in self.instructions:
            key = inst.name.lower()
            if key in {"barrier"}:
                continue
            counts[key] = counts.get(key, 0) + 1
        return counts

    def two_qubit_gate_count(self) -> int:
        return sum(1 for inst in self.instructions if inst.is_two_qubit())

    def depth(self) -> int:
        """Circuit depth using ASAP scheduling on logical qubits."""
        if not self.instructions:
            return 0

        qubit_layers: Dict[Qubit, int] = {}
        max_layer = 0
        for inst in self.instructions:
            if inst.name.lower() == "barrier":
                continue
            layers = [qubit_layers.get(q, 0) for q in inst.qubits]
            layer = max(layers, default=0)
            for q in inst.qubits:
                qubit_layers[q] = layer + 1
            max_layer = max(max_layer, layer + 1)
        return max_layer
