"""Reproduce the comparison of Escofet et al. (ACM TQC 2025) with qevedo-compiler.

    python benchmarks/multicore/run.py [--quick]

Three experiments, as in the paper (Section 6), each on QFT, the Draper and
Cuccaro adders, random circuits and Quantum Volume:

* virtual scaling: 10 cores of 10 qubits, circuits of 50 to 100 qubits;
* weak scaling: 200 qubits on 2 to 10 cores (as many qubits as slots);
* strong scaling: cores of 10 qubits, 2 to 20 cores, full.

and one the papers do not have:

* topology: 10 cores of 10 qubits, full, linked all to all, in a ring, in a
  2 × 5 grid or in a line, where moving between unlinked cores takes one
  teleportation per link. HQA counts links in its cost; "hqa-uniform-cost"
  is HQA with the papers' cost (every move is one), to show what that gains.

Every mapper starts from the same OEE placement. Results go to
results/multicore.csv; plot.py draws them.
"""

from __future__ import annotations

import argparse
import csv
import platform
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from circuits import BENCHMARKS, two_qubit_gates  # noqa: E402
from qevedo.compiler.multicore import map_to_cores  # noqa: E402

METHODS = ("hqa", "fgp-roee", "naive")
TOPOLOGY_METHODS = ("hqa", "hqa-uniform-cost", "fgp-roee", "naive")
OUT = Path(__file__).parent / "results" / "multicore.csv"


def experiments(quick: bool):
    """(experiment, cores, qubits per core, circuit qubits, topology) to run."""
    virtual = [("virtual", 10, 10, q, "all-to-all") for q in range(50, 101, 10)]
    weak = [("weak", c, 200 // c, 200, "all-to-all") for c in (2, 4, 5, 8, 10)]
    strong = [("strong", c, 10, 10 * c, "all-to-all") for c in range(2, 21, 2)]
    topology = [("topology", 10, 10, 100, t) for t in ("all-to-all", "ring", "grid:2x5", "line")]
    if quick:
        return virtual[::3] + weak[::2] + strong[::4] + topology
    return virtual + weak + strong + topology


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--quick", action="store_true", help="a few configurations only")
    args = parser.parse_args()
    rows = []
    for experiment, cores, per_core, qubits, topology in experiments(args.quick):
        methods = TOPOLOGY_METHODS if experiment == "topology" else METHODS
        for benchmark in BENCHMARKS:
            n, pairs = two_qubit_gates(benchmark, qubits)
            if n > cores * per_core or n < 2:
                continue
            for method in methods:
                mapping = map_to_cores(
                    pairs, cores, per_core, method=method, num_qubits=n, topology=topology
                )
                rows.append(
                    {
                        "experiment": experiment,
                        "benchmark": benchmark,
                        "cores": cores,
                        "qubits_per_core": per_core,
                        "topology": topology,
                        "qubits": n,
                        "two_qubit_gates": len(pairs),
                        "timeslices": len(mapping.timeslices),
                        "method": method,
                        "communications": mapping.communications,
                        "moved_qubits": mapping.moved_qubits,
                        "seconds": round(mapping.seconds, 6),
                    }
                )
                print(
                    f"{experiment:8} {topology:10} {benchmark:15} {cores:3}×{per_core:<3} {n:4} qubits "
                    f"{method:16} {mapping.communications:7} comms {mapping.seconds:8.3f} s",
                    flush=True,
                )
    OUT.parent.mkdir(exist_ok=True)
    with OUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nwrote {OUT} ({len(rows)} rows) on {platform.processor() or platform.machine()}")


if __name__ == "__main__":
    main()
