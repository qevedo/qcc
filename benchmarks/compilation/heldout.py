"""Held-out Hamiltonians on the chip: random Hamiltonians not used for tuning.

    pip install qiskit       # besides qevedo-compiler
    python benchmarks/compilation/heldout.py [--basis cx|cz] [--no-verify]

The chip Pauli-network changes (0.4.15 onwards) were tuned on the two
benchmark Hamiltonians of circuits.py (seeds 2008 and 2016) and judged on
these: one Trotter step of a random Hamiltonian built as
circuits.hamiltonian, n = 6..16 qubits with seeds 100n + 7 and 100n + 8
(22 circuits), compiled by QCC alone onto device.py's 19-qubit heavy-hex
chip. It prints two-qubit gates, gates in all and seconds per circuit and
the totals; each output is checked against the input's state as in
device.py unless --no-verify.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
from qiskit import QuantumCircuit, qasm2, transpile
from qiskit.circuit.library import PauliEvolutionGate
from qiskit.quantum_info import SparsePauliOp
from qiskit.transpiler import CouplingMap

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))

from circuits import INPUT_GATES  # noqa: E402
from device import BASES, qcc, same_state  # noqa: E402


def hamiltonian(n: int, seed: int, terms_per_qubit: int = 4) -> tuple[QuantumCircuit, str]:
    """circuits.hamiltonian with another seed, in the standard gates and as
    OpenQASM 2."""
    rng = np.random.default_rng(seed)
    labels = []
    for _ in range(terms_per_qubit * n):
        qubits = rng.choice(n, rng.integers(2, min(n, 6) + 1), replace=False)
        label = ["I"] * n
        for q in qubits:
            label[q] = "XYZ"[rng.integers(3)]
        labels.append("".join(label))
    op = SparsePauliOp(labels, coeffs=rng.uniform(-1, 1, len(labels)))
    qc = QuantumCircuit(n)
    qc.append(PauliEvolutionGate(op, time=0.3), range(n))
    qc = transpile(qc, basis_gates=INPUT_GATES, optimization_level=0)
    return qc, qasm2.dumps(qc)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--basis", choices=sorted(BASES), default="cx")
    parser.add_argument("--no-verify", action="store_true")
    args = parser.parse_args()
    coupling = CouplingMap.from_heavy_hex(3)
    two = total = 0
    seconds = 0.0
    failed = 0
    for n in range(6, 17):
        for seed in (100 * n + 7, 100 * n + 8):
            qc, qasm = hamiltonian(n, seed)
            out, final, time = qcc(qasm, BASES[args.basis], coupling)
            ops = [i.operation for i in out.data if i.operation.name not in ("barrier", "measure")]
            t2 = sum(1 for op in ops if op.num_qubits == 2)
            ok = True if args.no_verify else same_state(qc, out, final)
            failed += not ok
            two, total, seconds = two + t2, total + len(ops), seconds + time
            print(f"{args.basis} n={n:2} seed={seed} {t2:5} 2q {len(ops):6} gates {time:6.2f} s {ok}", flush=True)
    print(f"total {args.basis}: {two} two-qubit gates, {total} gates, {seconds:.1f} s; {failed} failed the check")
    if failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
