"""Qiskit with and without QCC as its optimization stage, on a routed chip.

    pip install qiskit       # besides qevedo-compiler
    python benchmarks/compilation/plugin.py

The circuits of circuits.py (8 and 16 qubits) are transpiled by Qiskit at
optimization level 3 onto a 19-qubit heavy-hex chip, for {rz, sx, x, cx} and
{rz, sx, x, cz}: once as is, once with optimization_method="qcc", which keeps
Qiskit's layout, routing and translation and replaces its optimization stage
with QCC's. Every output is checked: from |0…0⟩, it must reach the input's
state, on the physical qubits where the layout leaves each circuit qubit.
Results go to results/plugin.csv.
"""

from __future__ import annotations

import csv
import statistics
import sys
import time
from pathlib import Path

import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit.quantum_info import Statevector
from qiskit.transpiler import CouplingMap

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))

from circuits import FAMILIES, build  # noqa: E402

BASES = {"cx": ["rz", "sx", "x", "cx"], "cz": ["rz", "sx", "x", "cz"]}


def same_state(qc: QuantumCircuit, out: QuantumCircuit) -> bool:
    """Whether ``out`` maps |0…0⟩ to ``qc``'s state, up to a global phase, with
    circuit qubit v ending on physical qubit ``final_index_layout()[v]``."""
    layout = out.layout.final_index_layout()
    psi = Statevector(qc).data
    expected = np.zeros(2**out.num_qubits, dtype=complex)
    for index, amplitude in enumerate(psi):
        physical = sum(1 << layout[v] for v in range(qc.num_qubits) if index >> v & 1)
        expected[physical] = amplitude
    return abs(np.vdot(expected, Statevector(out).data)) > 1 - 1e-8


def main() -> None:
    coupling = CouplingMap.from_heavy_hex(3)
    rows = []
    for basis_name, basis in BASES.items():
        for family in FAMILIES:
            for n in (8, 16):
                qc, _ = build(family, n)
                for method in ("default", "qcc"):
                    start = time.perf_counter()
                    out = transpile(
                        qc,
                        coupling_map=coupling,
                        basis_gates=basis,
                        optimization_level=3,
                        optimization_method=None if method == "default" else method,
                        seed_transpiler=1,
                    )
                    seconds = time.perf_counter() - start
                    ops = [i.operation for i in out.data if i.operation.name not in ("barrier", "measure")]
                    verified = same_state(qc, out)
                    rows.append(
                        {
                            "basis": basis_name,
                            "family": family,
                            "qubits": qc.num_qubits,
                            "optimization": method,
                            "two_qubit_gates": sum(1 for op in ops if op.num_qubits == 2),
                            "gates": len(ops),
                            "seconds": round(seconds, 4),
                            "verified": verified,
                        }
                    )
                    r = rows[-1]
                    print(
                        f"{basis_name} {family:15} {r['qubits']:3} {method:8} {r['two_qubit_gates']:6} 2q "
                        f"{r['gates']:7} gates {seconds:7.3f} s {verified}",
                        flush=True,
                    )
    out = HERE / "results" / "plugin.csv"
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    by = {(r["basis"], r["family"], r["qubits"], r["optimization"]): r for r in rows}
    points = sorted({k[:3] for k in by})
    for key in ("two_qubit_gates", "gates"):
        ratio = statistics.geometric_mean(
            max(by[p + ("qcc",)][key], 1) / max(by[p + ("default",)][key], 1) for p in points
        )
        print(f"{key}: Qiskit + QCC / Qiskit = {ratio:.3f}x")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
