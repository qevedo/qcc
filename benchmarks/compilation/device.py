"""Compiling for a chip: QCC alone, Qiskit, and Qiskit with QCC's optimizations.

    pip install qiskit       # besides qevedo-compiler
    python benchmarks/compilation/device.py

The circuits of circuits.py (8 and 16 qubits) are compiled onto a 19-qubit
heavy-hex chip, for {rz, sx, x, cx} and {rz, sx, x, cz}:

- qcc: QCC alone (compile_circuit with the chip's couplings): it optimizes
  on logical qubits, lays out and routes with its SABRE, and optimizes again;
- qiskit: transpile at optimization level 3;
- qiskit+qcc: transpile with optimization_method="qcc" (see plugin.py).

Every output is checked: from |0…0⟩, it must reach the input's state, on the
physical qubits where each circuit qubit ends. Results go to
results/device.csv.
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
METHODS = ("qcc", "qiskit", "qiskit+qcc")


def same_state(qc: QuantumCircuit, out: QuantumCircuit, final: list[int]) -> bool:
    """Whether ``out`` maps |0…0⟩ to ``qc``'s state, up to a global phase, with
    circuit qubit v ending on physical qubit ``final[v]``."""
    psi = Statevector(qc).data
    expected = np.zeros(2**out.num_qubits, dtype=complex)
    for index, amplitude in enumerate(psi):
        physical = sum(1 << final[v] for v in range(qc.num_qubits) if index >> v & 1)
        expected[physical] = amplitude
    return abs(np.vdot(expected, Statevector(out).data)) > 1 - 1e-8


def qcc(qasm: str, basis: list[str], coupling: CouplingMap):
    from qevedo.compiler import CompileOptions, compile_circuit, emit_qasm
    from qevedo.compiler.device import DeviceSpec
    from qevedo.compiler.frontend.qasm import parse_qasm

    edges = sorted({tuple(sorted(e)) for e in coupling.get_edges()})
    device = DeviceSpec(
        name="heavy-hex", topology="heavy_hex", qubits=coupling.size(), edges=edges,
        native_gates=set(basis),
    )
    circuit = parse_qasm(source=qasm)
    start = time.perf_counter()
    out = compile_circuit(circuit, CompileOptions(device=device))
    seconds = time.perf_counter() - start
    final = [out.output_permutation[q].index for q in sorted(out.output_permutation, key=lambda q: q.index)]
    return QuantumCircuit.from_qasm_str(emit_qasm(out)), final, seconds


def qiskit(qc: QuantumCircuit, basis: list[str], coupling: CouplingMap, method: str | None):
    start = time.perf_counter()
    out = transpile(
        qc, coupling_map=coupling, basis_gates=basis, optimization_level=3,
        optimization_method=method, seed_transpiler=1,
    )
    seconds = time.perf_counter() - start
    return out, out.layout.final_index_layout(), seconds


def main() -> None:
    coupling = CouplingMap.from_heavy_hex(3)
    rows = []
    for basis_name, basis in BASES.items():
        for family in FAMILIES:
            for n in (8, 16):
                qc, qasm = build(family, n)
                for method in METHODS:
                    if method == "qcc":
                        out, final, seconds = qcc(qasm, basis, coupling)
                    else:
                        plugin = "qcc" if method == "qiskit+qcc" else None
                        out, final, seconds = qiskit(qc, basis, coupling, plugin)
                    ops = [i.operation for i in out.data if i.operation.name not in ("barrier", "measure")]
                    row = {
                        "basis": basis_name,
                        "family": family,
                        "qubits": qc.num_qubits,
                        "compiler": method,
                        "two_qubit_gates": sum(1 for op in ops if op.num_qubits == 2),
                        "gates": len(ops),
                        "seconds": round(seconds, 4),
                        "verified": same_state(qc, out, final),
                    }
                    rows.append(row)
                    print(
                        f"{basis_name} {family:15} {n:3} {method:11} {row['two_qubit_gates']:6} 2q "
                        f"{row['gates']:7} gates {seconds:7.3f} s {row['verified']}",
                        flush=True,
                    )
    out = HERE / "results" / "device.csv"
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    by = {(r["basis"], r["family"], r["qubits"], r["compiler"]): r for r in rows}
    print("\n| gates | compiler | two-qubit gates / Qiskit | gates / Qiskit | median time |")
    print("|---|---|---|---|---|")
    for basis_name in BASES:
        points = sorted({k[1:3] for k in by if k[0] == basis_name})
        for method in METHODS:
            ratio = lambda key: statistics.geometric_mean(  # noqa: E731
                max(by[(basis_name, *p, method)][key], 1) / max(by[(basis_name, *p, "qiskit")][key], 1)
                for p in points
            )
            time_ = statistics.median(by[(basis_name, *p, method)]["seconds"] for p in points)
            print(
                f"| {basis_name} | {method} | {ratio('two_qubit_gates'):.3f}× | {ratio('gates'):.3f}× "
                f"| {time_ * 1000:.0f} ms |"
            )
    wrong = [r for r in rows if not r["verified"]]
    print(f"\n{len(wrong)} outputs failed the check" + "".join(f"\n  {r}" for r in wrong))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
