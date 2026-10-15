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

and one more that takes the circuits down to physical qubits:

* routing: 100 qubits on 10 cores of 10, with the qubits of each core coupled
  all to all, in a 2 × 5 grid or in a line. Each core is routed with SABRE
  between teleportations (map_and_route), counting SWAPs as well.
* latency: the same circuits, cores in a 2 × 5 grid, linked all to all or in
  a line, scheduled with 1, 2, 4 or unlimited channels (EPR pairs at once)
  per link, with the default Timing otherwise; "hqa-latency" (HQA weighing
  how much each move delays its gate) is mapped for each capacity.
* ablation: the routing experiment's grid and line chips without choosing
  the slot of an arriving qubit by its next partner, and with SABRE looking
  0 to 64 timeslices past each epoch.

Every mapper starts from the same OEE placement. Results go to
results/{multicore,routing,latency,ablation}.csv; plot.py draws the
first.
"""

from __future__ import annotations

import argparse
import csv
import platform
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from circuits import BENCHMARKS, two_qubit_gates  # noqa: E402
from qevedo.compiler.multicore import Timing, map_and_route, map_to_cores  # noqa: E402

METHODS = ("hqa", "fgp-roee", "naive")
TOPOLOGY_METHODS = ("hqa", "hqa-uniform-cost", "fgp-roee", "naive")
LATENCY_METHODS = ("hqa", "hqa-latency", "fgp-roee", "naive")
OUT = Path(__file__).parent / "results" / "multicore.csv"
ROUTING_OUT = Path(__file__).parent / "results" / "routing.csv"
LATENCY_OUT = Path(__file__).parent / "results" / "latency.csv"
ABLATION_OUT = Path(__file__).parent / "results" / "ablation.csv"
#: (slot by next partner, timeslices SABRE looks past an epoch); (True, 4) is the default.
ABLATIONS = ((False, 0), (True, 0), (True, 1), (True, 2), (True, 4), (True, 8), (True, 16), (True, 64))
#: Channels per link; None is unlimited.
CAPACITIES = (1, 2, 4, None)
CHIPS = ("all-to-all", "grid:2x5", "line")


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
    write(OUT, rows)
    write(ROUTING_OUT, routing())
    write(LATENCY_OUT, latency())
    write(ABLATION_OUT, ablation())


def routing() -> list[dict]:
    """The routing experiment: teleportations and SWAPs per method and chip."""
    rows = []
    for chip in CHIPS:
        for benchmark in BENCHMARKS:
            n, pairs = two_qubit_gates(benchmark, 100)
            for method in METHODS:
                result = map_and_route(pairs, 10, 10, chip=chip, method=method, num_qubits=n)
                rows.append(
                    {
                        "benchmark": benchmark,
                        "chip": chip,
                        "qubits": n,
                        "two_qubit_gates": len(pairs),
                        "method": method,
                        "teleportations": result.teleportations,
                        "swaps": result.swaps,
                        "seconds": round(result.seconds, 6),
                    }
                )
                print(
                    f"routing  {chip:10} {benchmark:15} {method:16} {result.teleportations:7} "
                    f"teleportations {result.swaps:7} swaps {result.seconds:8.3f} s",
                    flush=True,
                )
    return rows


def latency() -> list[dict]:
    """The latency experiment: how long each mapping takes per link capacity."""
    rows = []
    for topology in ("all-to-all", "line"):
        for benchmark in BENCHMARKS:
            n, pairs = two_qubit_gates(benchmark, 100)
            for method in LATENCY_METHODS:
                result = None
                for capacity in CAPACITIES:
                    timing = Timing(link_capacity=capacity or len(pairs))
                    # hqa-latency plans for the timing, so it maps again per capacity.
                    if result is None or method == "hqa-latency":
                        result = map_and_route(
                            pairs,
                            10,
                            10,
                            chip="grid:2x5",
                            method=method,
                            num_qubits=n,
                            topology=topology,
                            timing=timing,
                        )
                    schedule = result.schedule(timing)
                    rows.append(
                        {
                            "benchmark": benchmark,
                            "topology": topology,
                            "chip": "grid:2x5",
                            "method": method,
                            "link_capacity": capacity or "unlimited",
                            "teleportations": result.teleportations,
                            "swaps": result.swaps,
                            "duration": schedule.duration,
                            "epr_wait": schedule.epr_wait,
                        }
                    )
                print(
                    f"latency  {topology:10} {benchmark:15} {method:16} "
                    + " ".join(f"{r['duration']:8.0f}" for r in rows[-len(CAPACITIES) :]),
                    flush=True,
                )
    return rows


def ablation() -> list[dict]:
    """The ablation experiment: SWAPs without the routing choices."""
    rows = []
    for chip in CHIPS[1:]:
        for benchmark in BENCHMARKS:
            n, pairs = two_qubit_gates(benchmark, 100)
            for method in METHODS:
                for partner_slots, look_ahead in ABLATIONS:
                    result = map_and_route(
                        pairs,
                        10,
                        10,
                        chip=chip,
                        method=method,
                        num_qubits=n,
                        partner_slots=partner_slots,
                        look_ahead=look_ahead,
                    )
                    rows.append(
                        {
                            "benchmark": benchmark,
                            "chip": chip,
                            "method": method,
                            "partner_slots": partner_slots,
                            "look_ahead": look_ahead,
                            "swaps": result.swaps,
                        }
                    )
                print(
                    f"ablation {chip:10} {benchmark:15} {method:16} "
                    + " ".join(f"{r['swaps']:6}" for r in rows[-len(ABLATIONS) :]),
                    flush=True,
                )
    return rows


def write(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nwrote {path} ({len(rows)} rows) on {platform.processor() or platform.machine()}")


if __name__ == "__main__":
    main()
