"""EPR pairs and SWAPs of TeleSABRE and of qevedo's mappers on the same machine and circuits.

    git clone https://github.com/Haimrich/telesabre && cd telesabre
    gcc -O3 -flto -I src src/*.c -o telesabre -lm
    python benchmarks/telesabre/compare.py path/to/telesabre [--rerun]

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
attempt fails is recorded with success False.

Our mappers run twice. With abstract communication (the method's name), a
qubit teleports from any slot to any free one, so SWAPs are not comparable
with TeleSABRE's. With ports ("+ports"), each link meets each core at the
same physical qubit as in device H: a qubit is routed onto it to teleport
and lands on the other core's, and the mapper keeps a slot free in every
core, as TeleSABRE does; SWAPs are then comparable too (TeleSABRE also
needs its communication qubits free, and can use telegates, which we do
not). Results go to results/telesabre.csv; TeleSABRE's rows are reused
from it unless --rerun is given, since TeleSABRE takes a quarter of an hour.
"""

from __future__ import annotations

import csv
import json
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
OUT = HERE / "results" / "telesabre.csv"


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
    seconds = round(time.perf_counter() - start, 3)
    if "No successful runs" in out:  # every attempt deadlocked
        return {
            "epr_pairs": "",
            "teleportations": "",
            "telegates": "",
            "swaps": "",
            "success": False,
            "seconds": seconds,
        }
    result = out[out.rfind("Result:") :]
    value = {k: int(re.search(rf"{k}: (\d+)", result).group(1)) for k in ("Teledata", "Telegate", "Swaps")}
    return {
        "epr_pairs": value["Teledata"] + value["Telegate"],
        "teleportations": value["Teledata"],
        "telegates": value["Telegate"],
        "swaps": value["Swaps"],
        "success": "Success: true" in result,
        "seconds": seconds,
    }


def circuits(root: Path, scratch: Path) -> list[tuple[str, str, Path]]:
    """(suite, name, file): TeleSABRE's 64-qubit circuits, then ours as OpenQASM."""
    out = [
        ("telesabre", p.name.split("_")[0], p)
        for p in sorted((root / "circuits" / "qasm_64").glob("*.qasm"))
    ]
    for name in OURS:
        n, pairs = two_qubit_gates(name, 64)
        if n > 64:  # the adders need more qubits than their size
            n, pairs = two_qubit_gates(name, 62)
        path = scratch / f"{name}_{n}.qasm"
        path.write_text(
            f'OPENQASM 2.0;\ninclude "qelib1.inc";\nqreg q[{n}];\n'
            + "".join(f"cx q[{a}],q[{b}];\n" for a, b in pairs)
        )
        out.append(("qevedo", name, path))
    return out


def device_ports(root: Path) -> dict[tuple[int, int], tuple[int, int]]:
    """Device H's links as (core, core) -> (slot, slot); cores have 16 qubits."""
    device = json.loads((root / "devices" / "H_grid_2_3_4_4.json").read_text())["device"]
    return {(p // 16, q // 16): (p % 16, q % 16) for p, q in device["inter_core_edges"]}


def cached_telesabre(rerun: bool) -> dict[tuple[str, str], dict]:
    if rerun or not OUT.exists():
        return {}
    with OUT.open() as handle:
        return {
            (r["suite"], r["circuit"]): r for r in csv.DictReader(handle) if r["method"] == "telesabre"
        }


def main() -> None:
    args = [a for a in sys.argv[1:] if a != "--rerun"]
    if len(args) != 1:
        sys.exit(__doc__)
    root = Path(args[0]).resolve()
    cached = cached_telesabre("--rerun" in sys.argv)
    ports = device_ports(root)
    rows = []
    with tempfile.TemporaryDirectory() as scratch:
        for suite, name, qasm in circuits(root, Path(scratch)):
            row = cached.get((suite, name)) or {
                "circuit": name,
                "suite": suite,
                "method": "telesabre",
                **telesabre(root, qasm),
            }
            rows.append(row)
            for port_map in (None, ports):
                for method in METHODS:
                    r = map_and_route(
                        qasm.read_text(),
                        cores=6,
                        qubits_per_core=16,
                        chip="grid:4x4",
                        topology="grid:2x3",
                        method=method,
                        ports=port_map,
                    )
                    rows.append(
                        {
                            "circuit": name,
                            "suite": suite,
                            "method": method + ("+ports" if port_map else ""),
                            "epr_pairs": r.teleportations,
                            "teleportations": r.teleportations,
                            "telegates": 0,
                            "swaps": r.swaps,
                            "success": True,
                            "seconds": round(r.seconds, 3),
                        }
                    )
            print(
                f"{name:15} "
                + "  ".join(
                    f"{r['method']}: {r['epr_pairs']}/{r['swaps']}" for r in rows[-1 - 2 * len(METHODS) :]
                ),
                flush=True,
            )
    OUT.parent.mkdir(exist_ok=True)
    with OUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nwrote {OUT} (EPR pairs/SWAPs per method)")


if __name__ == "__main__":
    main()
