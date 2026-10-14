"""Charts and a summary table from results/multicore.csv.

    python benchmarks/multicore/plot.py

Writes charts/{virtual,weak,strong}.svg (communications per benchmark, one
panel each) and prints the ratios against HQA as a Markdown table.
"""

from __future__ import annotations

import csv
import statistics
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).parent
STYLE = {
    "hqa": ("HQA (qevedo)", "#d1495b", "o"),
    "fgp-roee": ("FGP-rOEE", "#00798c", "s"),
    "naive": ("Naive", "#edae49", "^"),
}
TITLES = {
    "virtual": ("Virtual scaling: 10 cores × 10 qubits", "circuit qubits", "qubits"),
    "weak": ("Weak scaling: 200 qubits, full", "cores", "cores"),
    "strong": ("Strong scaling: cores of 10 qubits, full", "cores", "cores"),
}
NAMES = {
    "qft": "QFT",
    "draper": "Draper adder",
    "cuccaro": "Cuccaro adder",
    "random": "Random",
    "quantum_volume": "Quantum Volume",
}


def load():
    with (HERE / "results" / "multicore.csv").open() as handle:
        return list(csv.DictReader(handle))


def chart(rows, experiment):
    title, xlabel, xkey = TITLES[experiment]
    benchmarks = list(NAMES)
    fig, axes = plt.subplots(1, len(benchmarks), figsize=(17, 3.6))
    for ax, benchmark in zip(axes, benchmarks):
        for method, (label, color, marker) in STYLE.items():
            points = sorted(
                (int(r[xkey]), int(r["communications"]))
                for r in rows
                if r["experiment"] == experiment
                and r["benchmark"] == benchmark
                and r["method"] == method
            )
            if points:
                ax.plot(
                    *zip(*points),
                    marker=marker,
                    color=color,
                    label=label,
                    linewidth=1.5,
                    markersize=4,
                )
        ax.set_title(NAMES[benchmark], fontsize=10)
        ax.set_xlabel(xlabel)
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("non-local communications")
    axes[0].legend(fontsize=8)
    fig.suptitle(title + " (lower is better)", fontsize=11)
    fig.tight_layout()
    (HERE / "charts").mkdir(exist_ok=True)
    fig.savefig(HERE / "charts" / f"{experiment}.svg")
    plt.close(fig)


def summary(rows):
    by_point = defaultdict(dict)
    for r in rows:
        key = (r["experiment"], r["benchmark"], r["cores"], r["qubits"])
        by_point[key][r["method"]] = (int(r["communications"]), float(r["seconds"]))
    ratios = defaultdict(lambda: defaultdict(list))
    for (experiment, benchmark, _, _), m in by_point.items():
        hqa_comms, hqa_time = m["hqa"]
        for other in ("fgp-roee", "naive"):
            comms, _ = m[other]
            if hqa_comms:
                ratios[benchmark][other].append(comms / hqa_comms)
        ratios[benchmark]["time"].append(m["fgp-roee"][1] / max(hqa_time, 1e-6))
    print("| benchmark | FGP-rOEE / HQA communications | naive / HQA | FGP-rOEE / HQA time |")
    print("|---|---|---|---|")
    everything = defaultdict(list)
    for benchmark, name in NAMES.items():
        r = ratios[benchmark]
        for k in ("fgp-roee", "naive", "time"):
            everything[k] += r[k]
        print(
            f"| {name} | {statistics.median(r['fgp-roee']):.2f}× "
            f"| {statistics.median(r['naive']):.2f}× | {statistics.median(r['time']):.1f}× |"
        )
    # A few points where HQA takes almost no time inflate a mean of time
    # ratios, so the time column of the last row is a median.
    print(
        f"| **all points** | **{statistics.mean(everything['fgp-roee']):.2f}× (mean)** "
        f"| **{statistics.mean(everything['naive']):.2f}× (mean)** "
        f"| **{statistics.median(everything['time']):.1f}× (median)** |"
    )


if __name__ == "__main__":
    data = load()
    for experiment in TITLES:
        chart(data, experiment)
    summary(data)
