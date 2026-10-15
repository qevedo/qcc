# Paper draft

A preprint describing the multi-core toolkit: the reproduction of Escofet et
al.'s comparison, sparse links between cores, routing inside the cores, the
timing model and HQA-latency. Draft: the `[TODO: …]` marks in red are still
open (authors, the papers' own averages, timing values, further baselines,
an archived release).

```bash
python paper/make.py                 # tables, numbers and figures from the results
cd paper && latexmk -pdf main.tex    # or: tectonic main.tex
```

`make.py` reads `benchmarks/multicore/results/*.csv` and
`benchmarks/routing/results/sabre.csv` and writes `generated/`: a `.tex`
file per table, `numbers.tex` with every number the text quotes (as
macros) and the figures. Nothing in the paper is typed by hand from a
result; to refresh it, rerun the benchmarks
(`benchmarks/multicore/run.py`, `benchmarks/routing/vs_qiskit.py`) and
`make.py`.
