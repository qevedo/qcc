# Multi-core mapping benchmark

Reproduces the comparison of Escofet et al. (ACM TQC 2025): HQA against
FGP-rOEE and the naive mapper on QFT, the Draper and Cuccaro adders, random
circuits and Quantum Volume, in virtual, weak and strong scaling, and on
cores linked in a ring, a grid or a line (not in the papers).

```bash
pip install qiskit matplotlib     # besides qevedo-compiler
python benchmarks/multicore/run.py [--quick]
python benchmarks/multicore/plot.py
```

- `circuits.py` builds the circuits with Qiskit and keeps their CX gates
  (cached in `cache/`).
- `run.py` maps each circuit with each method and writes `results/multicore.csv`.
- `plot.py` draws `charts/` and prints the summary table.

Results and discussion: [docs/MULTICORE.md](../../docs/MULTICORE.md).
