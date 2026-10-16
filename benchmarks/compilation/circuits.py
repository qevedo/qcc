"""The circuits of the compilation benchmark, built with Qiskit.

Each family is a kind of circuit that compilers are commonly compared on;
each is written out as OpenQASM 2 with the standard gates of qelib1.inc, the
common input of every compiler. Everything is seeded, so the circuits are the
same on every run.
"""

from __future__ import annotations

import numpy as np
from qiskit import QuantumCircuit, qasm2, transpile
from qiskit.circuit.library import CDKMRippleCarryAdder, quantum_volume
from qiskit.synthesis.qft import synth_qft_full

# The standard gates the input may use: OpenQASM 2's qelib1.inc.
INPUT_GATES = ["u", "u1", "u2", "u3", "x", "y", "z", "h", "s", "sdg", "t", "tdg", "sx",
               "rx", "ry", "rz", "cx", "cz", "cy", "ch", "swap", "ccx", "crz", "cu1", "rxx", "rzz"]


def qft(n: int) -> QuantumCircuit:
    """The quantum Fourier transform, with its final SWAPs."""
    return synth_qft_full(n)


def qaoa(n: int, layers: int = 2) -> QuantumCircuit:
    """QAOA for MaxCut on a random 3-regular graph (n even)."""
    rng = np.random.default_rng(n)
    while True:  # pair up three copies of each vertex until no loops or repeats
        stubs = rng.permutation(np.repeat(np.arange(n), 3))
        edges = {tuple(sorted(stubs[i : i + 2])) for i in range(0, 3 * n, 2)}
        if len(edges) == 3 * n // 2 and all(a != b for a, b in edges):
            break
    qc = QuantumCircuit(n)
    qc.h(range(n))
    for layer in range(layers):
        gamma, beta = 0.4 + 0.1 * layer, 0.7 - 0.1 * layer
        for a, b in sorted(edges):
            qc.rzz(2 * gamma, int(a), int(b))
        qc.rx(2 * beta, range(n))
    return qc


def volume(n: int) -> QuantumCircuit:
    """Quantum Volume: n layers of random SU(4) blocks on random pairs."""
    return quantum_volume(n, seed=n)


def ghz(n: int) -> QuantumCircuit:
    qc = QuantumCircuit(n)
    qc.h(0)
    for q in range(n - 1):
        qc.cx(q, q + 1)
    return qc


def heisenberg(n: int, steps: int = 4) -> QuantumCircuit:
    """Trotterized Heisenberg chain: XX, YY and ZZ couplings on a line."""
    qc = QuantumCircuit(n)
    dt = 0.1
    for _ in range(steps):
        for parity in (0, 1):
            for q in range(parity, n - 1, 2):
                qc.rxx(2 * dt, q, q + 1)
                qc.ryy(2 * dt, q, q + 1)
                qc.rzz(2 * dt, q, q + 1)
    return qc


def adder(n: int) -> QuantumCircuit:
    """Cuccaro's ripple-carry adder of two (n − 2) / 2-bit numbers."""
    return CDKMRippleCarryAdder((n - 2) // 2, kind="full").decompose()


def clifford_t(n: int, depth: int = 40) -> QuantumCircuit:
    """Random layers of H, S, T and CX."""
    rng = np.random.default_rng(1000 + n)
    qc = QuantumCircuit(n)
    for _ in range(depth):
        order = rng.permutation(n)
        for k in range(0, n - 1, 2):
            if rng.random() < 0.5:
                qc.cx(int(order[k]), int(order[k + 1]))
        for q in range(n):
            gate = rng.integers(4)
            if gate == 0:
                qc.h(q)
            elif gate == 1:
                qc.s(q)
            elif gate == 2:
                qc.t(q)
    return qc


FAMILIES = {
    "qft": qft,
    "qaoa": qaoa,
    "quantum_volume": volume,
    "ghz": ghz,
    "heisenberg": heisenberg,
    "adder": adder,
    "clifford_t": clifford_t,
}


def build(family: str, n: int) -> tuple[QuantumCircuit, str]:
    """The circuit, in the standard gates, and as OpenQASM 2. Quantum Volume's
    random SU(4) blocks are given as Qiskit synthesizes them, in u3 and cx."""
    basis = ["u3", "cx"] if family == "quantum_volume" else INPUT_GATES
    circuit = FAMILIES[family](n)
    # One register, q: compilers order several registers differently.
    flat = QuantumCircuit(circuit.num_qubits)
    flat.compose(circuit, inplace=True)
    qc = transpile(flat, basis_gates=basis, optimization_level=0)
    return qc, qasm2.dumps(qc)


if __name__ == "__main__":
    for family in FAMILIES:
        qc, _ = build(family, 8)
        print(f"{family:15} {qc.num_qubits:3} qubits {dict(qc.count_ops())}")
