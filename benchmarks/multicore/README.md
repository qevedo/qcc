# Multi-core mapping benchmark

Reproduces the comparison of Escofet et al. (ACM TQC 2025): HQA against
FGP-rOEE and the naive mapper on QFT, the Draper and Cuccaro adders, random
circuits and Quantum Volume, in virtual, weak and strong scaling, and on
cores linked in a ring, a grid or a line (not in the papers), and, routing
inside each core, the SWAPs a circuit then needs when the qubits of a core
are coupled in a grid or a line, and how long it takes when links make EPR
pairs a few at a time.

```bash
pip install qiskit matplotlib     # besides qevedo-compiler
python benchmarks/multicore/run.py [--quick]
python benchmarks/multicore/plot.py
```

- `circuits.py` builds the circuits with Qiskit and keeps their CX gates
  (cached in `cache/`).
- `run.py` maps each circuit with each method and writes `results/multicore.csv`,
  then routes them inside the cores (`results/routing.csv`, and the
  ablation `results/ablation.csv`), and schedules them per link capacity
  (`results/latency.csv`), through ports (`results/latency_ports.csv`) and
  with faster and slower links (`results/timing.csv`).
- `plot.py` draws `charts/` and prints the summary table.

Results and discussion: [docs/MULTICORE.md](../../docs/MULTICORE.md).
