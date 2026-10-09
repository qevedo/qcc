"""Technology-agnostic backend description."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple

import yaml


@dataclass
class DeviceSpec:
    """Describes qubit connectivity and native gates for a target."""

    name: str
    topology: str
    qubits: int
    edges: List[Tuple[int, int]] = field(default_factory=list)
    native_gates: Set[str] = field(default_factory=set)
    directed_edges: bool = False
    max_parallel_2q: int = 1
    gate_rules: Dict[str, dict] = field(default_factory=dict)

    @classmethod
    def from_dict(cls, data: dict) -> DeviceSpec:
        edges = [tuple(edge) for edge in data.get("edges", [])]
        topology = data.get("topology", "coupling_graph")
        qubits = int(data["qubits"])
        if topology == "all_to_all" and not edges:
            edges = [(i, j) for i in range(qubits) for j in range(i + 1, qubits)]
        native = {gate.lower() for gate in data.get("native_gates", [])}
        return cls(
            name=data["name"],
            topology=topology,
            qubits=qubits,
            edges=edges,
            native_gates=native,
            directed_edges=bool(data.get("directed_edges", False)),
            max_parallel_2q=int(data.get("max_parallel_2q", 1)),
            gate_rules=data.get("gate_rules", {}),
        )

    @classmethod
    def from_yaml(cls, path: str | Path) -> DeviceSpec:
        with open(path, encoding="utf-8") as handle:
            data = yaml.safe_load(handle)
        return cls.from_dict(data)

    def adjacency(self) -> Dict[int, Set[int]]:
        graph: Dict[int, Set[int]] = {q: set() for q in range(self.qubits)}
        for a, b in self.edges:
            graph.setdefault(a, set()).add(b)
            if not self.directed_edges:
                graph.setdefault(b, set()).add(a)
        return graph

    def supports_gate(self, gate: str, qubits: Iterable[int]) -> bool:
        gate = gate.lower()
        if gate not in self.native_gates:
            return False
        if gate not in self.gate_rules:
            return True
        rule = self.gate_rules[gate]
        placement = rule.get("placement") or rule.get("on")
        if placement == "edge" and len(list(qubits)) == 2:
            pair = tuple(qubits)
            if self.directed_edges:
                return pair in self.edges
            return pair in self.edges or (pair[1], pair[0]) in self.edges
        return True


def default_device() -> DeviceSpec:
    """All-to-all superconducting-style device with a common IBM-like gateset."""
    qubits = 32
    edges = [(i, j) for i in range(qubits) for j in range(i + 1, qubits)]
    return DeviceSpec(
        name="all_to_all_32",
        topology="all_to_all",
        qubits=qubits,
        edges=edges,
        native_gates={"rz", "sx", "x", "cx"},
    )
