"""Tables, numbers and figures of the paper, from the benchmark results.

    python paper/make.py          # then: cd paper && latexmk -pdf main.tex

Reads benchmarks/multicore/results/*.csv and benchmarks/routing/results/
sabre.csv (run benchmarks/multicore/run.py and benchmarks/routing/
vs_qiskit.py to refresh them) and writes paper/generated/: one .tex file
per table, numbers.tex (every number the text quotes, as macros) and the
figures as PDF. Nothing in the paper is typed by hand from a result.
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
RESULTS = HERE.parent / "benchmarks" / "multicore" / "results"
OUT = HERE / "generated"
BENCHMARKS = {
    "qft": "QFT",
    "draper": "Draper adder",
    "cuccaro": "Cuccaro adder",
    "random": "Random",
    "quantum_volume": "Quantum Volume",
}
STYLE = {
    "hqa": ("HQA", "#d1495b", "o"),
    "hqa-latency": ("HQA-latency", "#30638e", "v"),
    "hqa-uniform-cost": ("HQA, uniform cost", "#8f2d56", "D"),
    "fgp-roee": ("FGP-rOEE", "#00798c", "s"),
    "naive": ("Naive", "#edae49", "^"),
}
CAPACITIES = ("1", "2", "4", "unlimited")


def load(path: Path) -> list[dict]:
    with path.open() as handle:
        return list(csv.DictReader(handle))


class Numbers:
    """The numbers the text quotes, written as LaTeX macros."""

    def __init__(self) -> None:
        self.lines: list[str] = []

    def add(self, name: str, value: str) -> None:
        self.lines.append(f"\\newcommand{{\\{name}}}{{{value}}}")

    def write(self) -> None:
        (OUT / "numbers.tex").write_text("\n".join(self.lines) + "\n")


def times(x: float) -> str:
    return f"{x:.2f}$\\times$"


def percent(x: float) -> str:
    return f"{100 * x:.0f}\\%"


def table(name: str, header: str, columns: str, rows: list[str]) -> None:
    body = "\n".join(rows)
    (OUT / f"{name}.tex").write_text(
        f"\\begin{{tabular}}{{{columns}}}\n\\toprule\n{header} \\\\\n\\midrule\n{body}\n"
        "\\bottomrule\n\\end{tabular}\n"
    )


def reproduction(numbers: Numbers) -> None:
    """Escofet et al.'s comparison: FGP-rOEE and naive against HQA."""
    rows = [r for r in load(RESULTS / "multicore.csv") if r["experiment"] != "topology"]
    points = defaultdict(dict)
    for r in rows:
        key = (r["experiment"], r["benchmark"], r["cores"], r["qubits"])
        points[key][r["method"]] = (int(r["communications"]), float(r["seconds"]))
    ratios = defaultdict(lambda: defaultdict(list))
    for (_, benchmark, _, _), m in points.items():
        comms, seconds = m["hqa"]
        for other in ("fgp-roee", "naive"):
            ratios[benchmark][other].append(m[other][0] / comms)
        ratios[benchmark]["time"].append(m["fgp-roee"][1] / max(seconds, 1e-6))
    lines = []
    everything = defaultdict(list)
    for benchmark, name in BENCHMARKS.items():
        r = ratios[benchmark]
        for k in ("fgp-roee", "naive", "time"):
            everything[k] += r[k]
        lines.append(
            f"{name} & {times(statistics.median(r['fgp-roee']))} & "
            f"{times(statistics.median(r['naive']))} & {statistics.median(r['time']):.1f}$\\times$ \\\\"
        )
    fgp = statistics.mean(everything["fgp-roee"])
    naive = statistics.mean(everything["naive"])
    speed = statistics.median(everything["time"])
    lines.append("\\midrule")
    lines.append(
        f"All points & {times(fgp)} (mean) & {times(naive)} (mean) & {speed:.1f}$\\times$ (median) \\\\"
    )
    table(
        "table-reproduction",
        "Benchmark & FGP-rOEE / HQA & Naive / HQA & FGP-rOEE / HQA time",
        "lrrr",
        lines,
    )
    numbers.add("points", str(len(points)))
    numbers.add("fgpOverHqa", times(fgp))
    numbers.add("naiveOverHqa", times(naive))
    numbers.add("hqaSpeedup", f"{speed:.1f}$\\times$")
    numbers.add("cuccaroFgpOverHqa", times(statistics.median(ratios["cuccaro"]["fgp-roee"])))


def links(numbers: Numbers) -> None:
    """Sparse links: communications counted in link hops."""
    rows = [r for r in load(RESULTS / "multicore.csv") if r["experiment"] == "topology"]
    value = {(r["benchmark"], r["topology"], r["method"]): int(r["communications"]) for r in rows}
    names = {"all-to-all": "All to all", "ring": "Ring", "grid:2x5": "2 $\\times$ 5 grid", "line": "Line"}
    lines = []
    for topology, name in names.items():
        ratio = {
            m: statistics.mean(value[(b, topology, m)] / value[(b, topology, "hqa")] for b in BENCHMARKS)
            for m in ("hqa-uniform-cost", "fgp-roee", "naive")
        }
        lines.append(
            f"{name} & {times(ratio['hqa-uniform-cost'])} & {times(ratio['fgp-roee'])} & "
            f"{times(ratio['naive'])} \\\\"
        )
        key = topology.split(":")[0].replace("-", "")
        numbers.add(f"uniform{key.capitalize()}", times(ratio["hqa-uniform-cost"]))
        numbers.add(f"fgp{key.capitalize()}", times(ratio["fgp-roee"]))
    table(
        "table-links",
        "Links & Uniform cost / HQA & FGP-rOEE / HQA & Naive / HQA",
        "lrrr",
        lines,
    )


def routing(numbers: Numbers) -> None:
    """Teleportations and SWAPs once each core is routed on its chip."""
    rows = load(RESULTS / "routing.csv")
    names = {"grid:2x5": "2 $\\times$ 5 grid", "line": "Line"}
    lines = []
    for chip, name in names.items():
        for k, method in enumerate(("hqa", "fgp-roee", "naive")):
            tel = sum(int(r["teleportations"]) for r in rows if r["chip"] == chip and r["method"] == method)
            swaps = sum(int(r["swaps"]) for r in rows if r["chip"] == chip and r["method"] == method)
            lines.append(f"{name if k == 0 else ''} & {STYLE[method][0]} & {tel} & {swaps} \\\\")
        if chip != "line":
            lines.append("\\midrule")
    table("table-routing", "Chip & Method & Teleportations & SWAPs", "llrr", lines)
    # Break-even: how many SWAPs a teleportation must cost for HQA's total to be lowest.
    get = {(r["chip"], r["benchmark"], r["method"]): (int(r["teleportations"]), int(r["swaps"])) for r in rows}
    worst = 0.0
    for chip in names:
        for b in BENCHMARKS:
            th, sh = get[(chip, b, "hqa")]
            for m in ("fgp-roee", "naive"):
                t, s = get[(chip, b, m)]
                if t > th and s < sh:
                    worst = max(worst, (sh - s) / (t - th))
    numbers.add("breakEven", f"{worst:.1f}")

    ablation = load(RESULTS / "ablation.csv")
    total = defaultdict(int)
    for r in ablation:
        total[(r["partner_slots"], int(r["look_ahead"]))] += int(r["swaps"])
    default = total[("True", 4)]
    horizons = [total[("True", h)] for h in (1, 2, 4, 8, 16, 64)]
    numbers.add("slotSaving", percent(1 - total[("True", 0)] / total[("False", 0)]))
    numbers.add("lookAheadSaving", percent(1 - default / total[("True", 0)]))
    numbers.add("horizonSpread", percent(max(horizons) / min(horizons) - 1))


def latency(numbers: Numbers) -> None:
    """Durations per link capacity, against HQA's."""
    rows = load(RESULTS / "latency.csv")
    d = defaultdict(dict)
    for r in rows:
        d[(r["topology"], r["link_capacity"], r["benchmark"])][r["method"]] = (
            float(r["duration"]),
            int(r["teleportations"]),
        )
    names = {"all-to-all": "All to all", "line": "Line"}
    lines = []
    means = defaultdict(dict)
    for topology, name in names.items():
        for k, c in enumerate(CAPACITIES):
            keys = [key for key in d if key[0] == topology and key[1] == c]
            ratio = {
                m: statistics.mean(d[key][m][0] / d[key]["hqa"][0] for key in keys)
                for m in ("hqa-latency", "fgp-roee", "naive")
            }
            means[topology][c] = ratio
            lines.append(
                f"{name if k == 0 else ''} & {c} & \\textbf{{{times(ratio['hqa-latency'])}}} & "
                f"{times(ratio['fgp-roee'])} & {times(ratio['naive'])} \\\\"
            )
        if topology != "line":
            lines.append("\\midrule")
    table(
        "table-latency",
        "Links & Channels & HQA-latency / HQA & FGP-rOEE / HQA & Naive / HQA",
        "llrrr",
        lines,
    )
    numbers.add("fgpSlowerOneChannelAll", times(means["all-to-all"]["1"]["fgp-roee"]))
    numbers.add("fgpSlowerOneChannelLine", times(means["line"]["1"]["fgp-roee"]))
    numbers.add("fgpUnlimitedAll", times(means["all-to-all"]["unlimited"]["fgp-roee"]))
    numbers.add("naiveUnlimitedAll", times(means["all-to-all"]["unlimited"]["naive"]))
    lat = [means[t][c]["hqa-latency"] for t in names for c in CAPACITIES]
    numbers.add("latencyMeanBest", percent(1 - min(lat)))
    numbers.add("latencyMeanWorst", percent(1 - max(lat)))
    per_run = [v["hqa-latency"][0] / v["hqa"][0] for v in d.values()]
    numbers.add("latencyBestRun", percent(1 - min(per_run)))
    numbers.add("latencyWorstRun", f"{max(per_run):.2f}")
    wins = sum(1 for v in d.values() if v["hqa-latency"][0] <= min(x[0] for x in v.values()) + 1e-9)
    numbers.add("latencyWins", str(wins))
    numbers.add("latencyRuns", str(len(d)))
    extra = max(
        sum(d[k]["hqa-latency"][1] for k in d if k[:2] == (t, c))
        / sum(d[k]["hqa"][1] for k in d if k[:2] == (t, c))
        for t in names
        for c in CAPACITIES
    )
    numbers.add("latencyExtraTeleportations", percent(extra - 1))
    latency_figure(means)


def latency_figure(means: dict) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(7, 2.6), sharey=True)
    for ax, (topology, title) in zip(axes, (("all-to-all", "Cores linked all to all"), ("line", "Cores in a line"))):
        xs = range(len(CAPACITIES))
        ax.axhline(1.0, color=STYLE["hqa"][1], linewidth=1.2, label="HQA")
        for method in ("hqa-latency", "fgp-roee", "naive"):
            label, color, marker = STYLE[method]
            ax.plot(xs, [means[topology][c][method] for c in CAPACITIES], marker=marker, color=color, label=label)
        ax.set_xticks(list(xs), ["1", "2", "4", "$\\infty$"])
        ax.set_xlabel("channels per link")
        ax.set_title(title, fontsize=9)
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("duration / HQA's")
    axes[1].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(OUT / "fig-latency.pdf", metadata={"CreationDate": None})
    plt.close(fig)


def scaling_figure() -> None:
    """Strong scaling: communications per benchmark."""
    rows = load(RESULTS / "multicore.csv")
    fig, axes = plt.subplots(1, len(BENCHMARKS), figsize=(7, 2.2))
    for ax, (benchmark, name) in zip(axes, BENCHMARKS.items()):
        for method in ("hqa", "fgp-roee", "naive"):
            label, color, marker = STYLE[method]
            points = sorted(
                (int(r["cores"]), int(r["communications"]))
                for r in rows
                if r["experiment"] == "strong" and r["benchmark"] == benchmark and r["method"] == method
            )
            ax.plot(*zip(*points), marker=marker, color=color, label=label, linewidth=1, markersize=2.5)
        ax.set_title(name, fontsize=7)
        ax.tick_params(labelsize=6)
        ax.set_xlabel("cores", fontsize=7)
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("communications", fontsize=7)
    axes[0].legend(fontsize=5)
    fig.tight_layout()
    fig.savefig(OUT / "fig-strong.pdf", metadata={"CreationDate": None})
    plt.close(fig)


def sabre(numbers: Numbers) -> None:
    rows = load(HERE.parent / "benchmarks" / "routing" / "results" / "sabre.csv")
    numbers.add("sabreOurs", str(sum(int(r["qevedo_swaps"]) for r in rows)))
    numbers.add("sabreQiskit", str(sum(int(r["qiskit_swaps"]) for r in rows)))
    numbers.add("sabreCircuits", str(len(rows)))


def telesabre(numbers: Numbers) -> None:
    """EPR pairs against TeleSABRE's reference implementation, on its device H."""
    rows = load(HERE.parent / "benchmarks" / "telesabre" / "results" / "telesabre.csv")
    names = {
        "ae": "Amplitude estimation",
        "ghz": "GHZ",
        "graphstate": "Graph state",
        "qnn": "QNN",
        **BENCHMARKS,
    }
    methods = ("telesabre", "hqa", "hqa-latency", "fgp-roee")
    by = defaultdict(dict)
    order = []
    for r in rows:
        key = (r["suite"], r["circuit"])
        if key not in by:
            order.append(key)
        by[key][r["method"]] = r
    lines = []
    ratios = []
    failed = 0
    for k, key in enumerate(order):
        suite, circuit = key
        if k and suite != order[k - 1][0]:
            lines.append("\\midrule")
        cells = []
        for m in methods:
            r = by[key][m]
            cells.append(r["epr_pairs"] if r["success"] == "True" else "fails")
        ts = by[key]["telesabre"]
        if ts["success"] == "True":
            ratios.append(int(by[key]["hqa"]["epr_pairs"]) / int(ts["epr_pairs"]))
        else:
            failed += 1
        source = "MQT Bench" if suite == "telesabre" else "ours"
        lines.append(f"{names[circuit]} ({source}) & " + " & ".join(cells) + " \\\\")
    table(
        "table-telesabre",
        "Circuit & TeleSABRE & HQA & HQA-latency & FGP-rOEE",
        "lrrrr",
        lines,
    )
    numbers.add("telesabreHqaRatio", times(statistics.geometric_mean(ratios)))
    numbers.add("telesabreFailed", str(failed))
    numbers.add("telesabreCircuits", str(len(order)))


def main() -> None:
    OUT.mkdir(exist_ok=True)
    numbers = Numbers()
    reproduction(numbers)
    links(numbers)
    routing(numbers)
    latency(numbers)
    sabre(numbers)
    telesabre(numbers)
    scaling_figure()
    numbers.write()
    print(f"wrote {OUT}: " + ", ".join(sorted(p.name for p in OUT.iterdir())))


if __name__ == "__main__":
    main()
