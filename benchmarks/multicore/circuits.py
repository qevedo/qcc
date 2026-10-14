"""The benchmark circuits of Escofet et al., built with Qiskit, as two-qubit gate lists.

Each circuit is transpiled to {cx, u} so every two-qubit gate is a CX, and
only the qubit pairs of the CX gates are kept: single-qubit gates do not
affect where qubits must be. Results are cached in ``cache/`` (not committed).
"""

from __future__ import annotations

import json
from pathlib import Path

from qiskit import QuantumCircuit, transpile
from qiskit.circuit.library import (
    CDKMRippleCarryAdder,
    DraperQFTAdder,
    QuantumVolume,
)
from qiskit.circuit.random import random_circuit
from qiskit.synthesis import synth_qft_full

CACHE = Path(__file__).parent / "cache"
BENCHMARKS = ("qft", "draper", "cuccaro", "random", "quantum_volume")


def build(name: str, qubits: int, seed: int = 1) -> QuantumCircuit:
    """Benchmark ``name`` on (at most) ``qubits`` qubits."""
    if name == "qft":
        return synth_qft_full(qubits, do_swaps=False)
    if name == "draper":
        # Two registers of qubits // 2.
        return DraperQFTAdder(qubits // 2, kind="fixed").decompose()
    if name == "cuccaro":
        # Two registers of k, carry in and carry out: 2k + 2 qubits.
        return CDKMRippleCarryAdder((qubits - 2) // 2, kind="full").decompose()
    if name == "random":
        # Depth twice the number of qubits, gates on one or two qubits.
        return random_circuit(qubits, 2 * qubits, max_operands=2, seed=seed)
    if name == "quantum_volume":
        return QuantumVolume(qubits, seed=seed).decompose()
    raise ValueError(f"unknown benchmark {name!r}")


def two_qubit_gates(name: str, qubits: int, seed: int = 1) -> tuple[int, list[tuple[int, int]]]:
    """The circuit's number of qubits and its CX gates as qubit pairs, cached."""
    path = CACHE / f"{name}_{qubits}_{seed}.json"
    if path.exists():
        data = json.loads(path.read_text())
        return data["qubits"], [tuple(p) for p in data["pairs"]]
    circuit = transpile(build(name, qubits, seed), basis_gates=["cx", "u"], optimization_level=0)
    pairs = [
        (circuit.find_bit(i.qubits[0]).index, circuit.find_bit(i.qubits[1]).index)
        for i in circuit.data
        if i.operation.name == "cx"
    ]
    CACHE.mkdir(exist_ok=True)
    path.write_text(json.dumps({"qubits": circuit.num_qubits, "pairs": pairs}))
    return circuit.num_qubits, pairs
