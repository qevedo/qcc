"""EPR pairs used by TeleSABRE and by qevedo's mappers on the same machine and circuits.

    git clone https://github.com/Haimrich/telesabre && cd telesabre
    gcc -O3 -flto -I src src/*.c -o telesabre -lm
    python benchmarks/telesabre/compare.py path/to/telesabre

TeleSABRE (Russo et al., IEEE QCE 2025) routes onto multi-core chips with
teleportations and remote gates (telegates) through communication qubits at
fixed places; each uses one EPR pair. Its device H is 6 cores of 4 × 4
qubits in a 2 × 3 grid of cores, one link between neighbouring cores. We
map the same circuits onto the same machine (topology "grid:2x3", chip
"grid:4x4") with map_and_route, whose teleportations count one EPR pair per
link crossed.

The circuits are TeleSABRE's own 64-qubit MQT Bench circuits and the
64-qubit circuits of benchmarks/multicore (CX gates only). TeleSABRE retries
with new seeds up to 100 times when it deadlocks; a circuit where every
attempt fails is recorded with success False. Only EPR pairs are compared: our model leaves communication qubits abstract (a qubit can
teleport from any slot), so it needs no SWAPs to reach them, and its SWAP
counts are not comparable with TeleSABRE's. Results go to
results/telesabre.csv.
"""

from __future__ import annotations

import csv
import re
import subprocess
import sys
import tempfile
import time
from pathlib import Path

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE.parent / "multicore"))

from circuits import two_qubit_gates  # noqa: E402
from qevedo.compiler.multicore import map_and_route  # noqa: E402

METHODS = ("hqa", "hqa-latency", "fgp-roee")
OURS = ("qft", "draper", "cuccaro", "random", "quantum_volume")


def telesabre(root: Path, qasm: Path) -> dict:
    """TeleSABRE's result on device H, with its default configuration."""
    start = time.perf_counter()
    out = subprocess.run(
        [str(root / "telesabre"), "configs/default.json", "devices/H_grid_2_3_4_4.json", str(qasm)],
        cwd=root,
        capture_output=True,
        text=True,
        timeout=3600,
    ).stdout
    seconds = time.perf_counter() - start
    if "No successful runs" in out:  # every attempt deadlocked
        return {
            "epr_pairs": "",
            "teleportations": "",
            "telegates": "",
            "swaps": "",
            "success": False,
            "seconds": round(seconds, 3),
        }
    result = out[out.rfind("Result:") :]
    value = {k: int(re.search(rf"{k}: (\d+)", result).group(1)) for k in ("Teledata", "Telegate", "Swaps")}
    return {
        "epr_pairs": value["Teledata"] + value["Telegate"],
        "teleportations": value["Teledata"],
        "telegates": value["Telegate"],
        "swaps": value["Swaps"],
        "success": "Success: true" in result,
        "seconds": round(seconds, 3),
    }


def circuits(root: Path, scratch: Path) -> list[tuple[str, Path]]:
    """TeleSABRE's 64-qubit circuits, then ours written as OpenQASM."""
    out = [(p.name.split("_")[0], p) for p in sorted((root / "circuits" / "qasm_64").glob("*.qasm"))]
    for name in OURS:
        n, pairs = two_qubit_gates(name, 64)
        if n > 64:  # the adders need more qubits than their size
            n, pairs = two_qubit_gates(name, 62)
        path = scratch / f"{name}_{n}.qasm"
        path.write_text(
            f'OPENQASM 2.0;\ninclude "qelib1.inc";\nqreg q[{n}];\n'
            + "".join(f"cx q[{a}],q[{b}];\n" for a, b in pairs)
        )
        out.append((name, path))
    return out


def main() -> None:
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    root = Path(sys.argv[1]).resolve()
    rows = []
    with tempfile.TemporaryDirectory() as scratch:
        for name, qasm in circuits(root, Path(scratch)):
            source = "telesabre" if "nativegates" in qasm.name else "qevedo"
            ts = telesabre(root, qasm)
            rows.append({"circuit": name, "suite": source, "method": "telesabre", **ts})
            for method in METHODS:
                r = map_and_route(
                    qasm.read_text(),
                    cores=6,
                    qubits_per_core=16,
                    chip="grid:4x4",
                    topology="grid:2x3",
                    method=method,
                )
                rows.append(
                    {
                        "circuit": name,
                        "suite": source,
                        "method": method,
                        "epr_pairs": r.teleportations,
                        "teleportations": r.teleportations,
                        "telegates": 0,
                        "swaps": r.swaps,
                        "success": True,
                        "seconds": round(r.seconds, 3),
                    }
                )
            print(
                f"{name:15} " + "  ".join(f"{r['method']}: {r['epr_pairs']}" for r in rows[-1 - len(METHODS) :]),
                flush=True,
            )
    path = HERE / "results" / "telesabre.csv"
    path.parent.mkdir(exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nwrote {path}")


if __name__ == "__main__":
    main()
