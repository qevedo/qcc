"""Compile the same circuits with QCC, Qiskit, pytket and UCC, and compare.

    pip install qiskit pytket ucc     # besides qevedo-compiler; UCC needs Python 3.12-3.13
    python benchmarks/compilation/run.py [--sizes 8 16 32]

Every compiler gets the same OpenQASM 2 input (circuits.py) and compiles it
for an IBM-style gate set, {rz, sx, x, cx}, with every qubit coupled to every
other (no routing):

- QCC: qevedo.compiler's default pipeline (the Rust core);
- Qiskit: transpile at optimization level 3;
- pytket: FullPeepholeOptimise, then a rebase to the gate set;
- UCC: ucc.compile, its default pipeline.

Qiskit and pytket may remove SWAP gates by relabelling the output qubits
(Qiskit's ElidePermutations, pytket's implicit wire swaps); QCC and UCC do
not. That is part of what each does by default, so it is kept, and the
checks below account for it.

Times exclude reading the input, for every compiler.

For circuits of up to 10 qubits, each output is checked against the input:
the unitaries must agree up to a global phase. Compilers that are not
installed are skipped. Results go to results/compilation.csv.
"""

from __future__ import annotations

import argparse
import csv
import statistics
import sys
import time
from pathlib import Path

from qiskit import QuantumCircuit, transpile
from qiskit.quantum_info import Operator

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))

from circuits import FAMILIES, build  # noqa: E402

GATES = ["rz", "sx", "x", "cx"]
CHECK_UP_TO = 10


def qcc(qasm: str, n: int):
    from qevedo.compiler import CompileOptions, compile_circuit, emit_qasm
    from qevedo.compiler.device import DeviceSpec
    from qevedo.compiler.frontend.qasm import parse_qasm

    device = DeviceSpec(name="bench", topology="all_to_all", qubits=n, native_gates=set(GATES))
    circuit = parse_qasm(source=qasm)
    start = time.perf_counter()
    out = compile_circuit(circuit, CompileOptions(device=device))
    seconds = time.perf_counter() - start
    counts = out.gate_counts()
    return out.two_qubit_gate_count(), sum(counts.values()), seconds, lambda: Operator(
        QuantumCircuit.from_qasm_str(emit_qasm(out))
    )


def qiskit(qasm: str, n: int):
    qc = QuantumCircuit.from_qasm_str(qasm)
    start = time.perf_counter()
    out = transpile(qc, basis_gates=GATES, optimization_level=3, seed_transpiler=1)
    seconds = time.perf_counter() - start
    return two_and_total(out) + (seconds, lambda: Operator.from_circuit(out))


def pytket(qasm: str, n: int):
    from pytket import OpType
    from pytket.passes import AutoRebase, FullPeepholeOptimise
    from pytket.qasm import circuit_from_qasm_str, circuit_to_qasm_str

    circuit = circuit_from_qasm_str(qasm)
    start = time.perf_counter()
    FullPeepholeOptimise().apply(circuit)
    AutoRebase({OpType.CX, OpType.Rz, OpType.SX, OpType.X}).apply(circuit)
    seconds = time.perf_counter() - start
    two, total = circuit.n_2qb_gates(), circuit.n_gates

    def unitary():
        circuit.replace_implicit_wire_swaps()
        return Operator(QuantumCircuit.from_qasm_str(circuit_to_qasm_str(circuit)))

    return two, total, seconds, unitary


def ucc(qasm: str, n: int):
    from ucc import compile as ucc_compile

    qc = QuantumCircuit.from_qasm_str(qasm)
    start = time.perf_counter()
    out = ucc_compile(qc, target_gateset=set(GATES))
    seconds = time.perf_counter() - start
    return two_and_total(out) + (seconds, lambda: Operator.from_circuit(out))


def two_and_total(qc: QuantumCircuit) -> tuple[int, int]:
    ops = [i.operation for i in qc.data if i.operation.name not in ("barrier", "measure")]
    return sum(1 for op in ops if op.num_qubits == 2), len(ops)


COMPILERS = {"qcc": qcc, "qiskit": qiskit, "pytket": pytket, "ucc": ucc}


def available() -> dict:
    out = {}
    for name, run in COMPILERS.items():
        try:
            run("OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[2];\ncx q[0],q[1];\n", 2)
            out[name] = run
        except ImportError:
            print(f"skipping {name}: not installed")
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--sizes", type=int, nargs="+", default=[8, 16, 32])
    args = parser.parse_args()
    compilers = available()
    rows = []
    for family in FAMILIES:
        for n in args.sizes:
            qc, qasm = build(family, n)
            reference = Operator(qc) if qc.num_qubits <= CHECK_UP_TO else None
            for name, run in compilers.items():
                two, total, seconds, unitary = run(qasm, qc.num_qubits)
                verified = "" if reference is None else reference.equiv(unitary())
                rows.append(
                    {
                        "family": family,
                        "qubits": qc.num_qubits,
                        "compiler": name,
                        "two_qubit_gates": two,
                        "gates": total,
                        "seconds": round(seconds, 4),
                        "verified": verified,
                    }
                )
                print(
                    f"{family:15} {qc.num_qubits:3} {name:7} {two:6} 2q {total:7} gates "
                    f"{seconds:8.3f} s {'' if verified == '' else 'ok' if verified else 'WRONG'}",
                    flush=True,
                )
    out = HERE / "results" / "compilation.csv"
    out.parent.mkdir(exist_ok=True)
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    summary(rows, list(compilers))
    print(f"\nwrote {out}")


def summary(rows: list[dict], compilers: list[str]) -> None:
    """Geometric means against Qiskit, and median times."""
    by = {(r["family"], r["qubits"], r["compiler"]): r for r in rows}
    points = sorted({(r["family"], r["qubits"]) for r in rows})
    print("\n| compiler | two-qubit gates / Qiskit | gates / Qiskit | median time |")
    print("|---|---|---|---|")
    for c in compilers:
        ratio = lambda key: statistics.geometric_mean(  # noqa: E731
            max(by[p + (c,)][key], 1) / max(by[p + ("qiskit",)][key], 1) for p in points
        )
        time_ = statistics.median(by[p + (c,)]["seconds"] for p in points)
        print(f"| {c} | {ratio('two_qubit_gates'):.3f}× | {ratio('gates'):.3f}× | {time_ * 1000:.0f} ms |")
    wrong = [r for r in rows if r["verified"] is False]
    print(f"\n{len(wrong)} outputs failed the check" + "".join(f"\n  {r}" for r in wrong))


if __name__ == "__main__":
    main()
