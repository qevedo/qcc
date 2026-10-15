"""SWAPs inserted by qevedo's SABRE router and by Qiskit's, on the same circuits and chips.

    pip install qiskit       # besides qevedo-compiler
    python benchmarks/routing/vs_qiskit.py

The circuits are the CX gates of the multi-core benchmarks (QFT, random,
Quantum Volume, Cuccaro adder), on a line of 20 qubits, a 5 × 5 grid and
a 19-qubit heavy-hex chip. Both routers search for the layout too: Qiskit's
SabreLayout with 8 layout trials, 8 swap trials and 4 iterations, ours with
its defaults (8 layout trials, 4 swap trials, 3 forward–backward passes).
Results go to results/sabre.csv.
"""

from __future__ import annotations

import csv
import sys
import time
from pathlib import Path

from qiskit import QuantumCircuit
from qiskit.transpiler import CouplingMap, PassManager
from qiskit.transpiler.passes import SabreLayout

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parent / "multicore"))

from circuits import two_qubit_gates  # noqa: E402
from qevedo.compiler import _core  # noqa: E402

BENCHMARKS = ("qft", "random", "quantum_volume", "cuccaro")


def grid(rows: int, columns: int) -> list[tuple[int, int]]:
    edges = []
    for r in range(rows):
        for c in range(columns):
            p = r * columns + c
            if c + 1 < columns:
                edges.append((p, p + 1))
            if r + 1 < rows:
                edges.append((p, p + columns))
    return edges


def chips() -> dict[str, tuple[int, list[tuple[int, int]]]]:
    heavy_hex = CouplingMap.from_heavy_hex(3)
    return {
        "line 20": (20, [(i, i + 1) for i in range(19)]),
        "grid 5x5": (25, grid(5, 5)),
        "heavy-hex 19": (heavy_hex.size(), sorted({tuple(sorted(e)) for e in heavy_hex.get_edges()})),
    }


def main() -> None:
    rows = []
    for chip, (n, edges) in chips().items():
        for benchmark in BENCHMARKS:
            qubits, pairs = two_qubit_gates(benchmark, n)
            if qubits > n:  # the adders need two qubits more than their size
                qubits, pairs = two_qubit_gates(benchmark, n - 2)
            circuit = QuantumCircuit(n)
            for a, b in pairs:
                circuit.cx(a, b)
            coupling = CouplingMap([list(e) for e in edges] + [[b, a] for a, b in edges])
            start = time.perf_counter()
            layout = SabreLayout(coupling, seed=1, max_iterations=4, swap_trials=8, layout_trials=8)
            routed = PassManager([layout]).run(circuit)
            qiskit_seconds = time.perf_counter() - start
            start = time.perf_counter()
            _, _, _, swaps = _core.route([("cx", [a, b], [], []) for a, b in pairs], n, n, edges)
            ours_seconds = time.perf_counter() - start
            rows.append(
                {
                    "chip": chip,
                    "benchmark": benchmark,
                    "qubits": qubits,
                    "two_qubit_gates": len(pairs),
                    "qiskit_swaps": routed.count_ops().get("swap", 0),
                    "qevedo_swaps": swaps,
                    "qiskit_seconds": round(qiskit_seconds, 4),
                    "qevedo_seconds": round(ours_seconds, 4),
                }
            )
            print(", ".join(f"{k}={v}" for k, v in rows[-1].items()), flush=True)
    out = HERE / "results" / "sabre.csv"
    out.parent.mkdir(exist_ok=True)
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    total = {k: sum(r[k] for r in rows) for k in ("qiskit_swaps", "qevedo_swaps")}
    print(f"\nwrote {out}; total SWAPs: Qiskit {total['qiskit_swaps']}, qevedo {total['qevedo_swaps']}")


if __name__ == "__main__":
    main()
