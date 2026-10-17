"""Check QCC against its recorded benchmark results, for CI.

    pip install qiskit      # besides qevedo-compiler
    python benchmarks/compilation/check.py >> "$GITHUB_STEP_SUMMARY"

Compiles every circuit of the compilation benchmark with QCC alone, checks
the outputs of up to 10 qubits against the input's unitary, and compares the
gate counts with QCC's rows in results/compilation.csv (written by run.py,
next to Qiskit's, pytket's and UCC's). Prints a Markdown report, and exits
with an error if an output is wrong or uses more two-qubit gates, or more
gates in all, than recorded. When QCC gets better, rerun run.py to record it.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

from qiskit.quantum_info import Operator

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))

from circuits import build  # noqa: E402
from run import CHECK_UP_TO, qcc  # noqa: E402


def main() -> None:
    with (HERE / "results" / "compilation.csv").open() as handle:
        recorded = list(csv.DictReader(handle))
    qiskit = {(r["family"], r["qubits"]): r for r in recorded if r["compiler"] == "qiskit"}
    failures = []
    lines = [
        "## QCC compilation benchmark",
        "",
        "| circuit | qubits | two-qubit gates | gates | recorded | Qiskit | time |",
        "|---|---|---|---|---|---|---|",
    ]
    for row in (r for r in recorded if r["compiler"] == "qcc"):
        family, n = row["family"], int(row["qubits"])
        qc, qasm = build(family, n)
        two, total, seconds, unitary = qcc(qasm, qc.num_qubits)
        if qc.num_qubits <= CHECK_UP_TO and not Operator(qc).equiv(unitary()):
            failures.append(f"{family} {n}: the output is not equivalent to the input")
        before = (int(row["two_qubit_gates"]), int(row["gates"]))
        if (two, total) > before or total > before[1]:
            failures.append(f"{family} {n}: {two} / {total} gates, recorded {before[0]} / {before[1]}")
        q = qiskit[(family, row["qubits"])]
        lines.append(
            f"| {family} | {n} | {two} | {total} | {before[0]} / {before[1]} "
            f"| {q['two_qubit_gates']} / {q['gates']} | {seconds * 1000:.0f} ms |"
        )
    lines.append("")
    lines += [f"- ❌ {f}" for f in failures] or ["All outputs checked and no worse than recorded."]
    print("\n".join(lines))
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
