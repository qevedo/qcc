# Compilation benchmark

QCC against Qiskit, pytket and UCC on the same circuits, compiled for an
IBM-style gate set, {rz, sx, x, cx}, with every qubit coupled to every other
(no routing).

```bash
pip install qiskit pytket ucc     # besides qevedo-compiler; UCC needs Python 3.12–3.13
python benchmarks/compilation/run.py [--sizes 8 16 32]
```

- `circuits.py` builds seven families with Qiskit, seeded: QFT, QAOA (MaxCut
  on a random 3-regular graph, two layers), Quantum Volume, GHZ, a
  Trotterized Heisenberg chain, Cuccaro's adder and random Clifford+T
  layers. Each goes to every compiler as the same OpenQASM 2.
- `run.py` compiles them with QCC's default pipeline (relabelling outputs,
  as Qiskit and pytket do), Qiskit's `transpile` at optimization level 3,
  pytket's `FullPeepholeOptimise` plus a rebase, and UCC's `compile`, and
  writes `results/compilation.csv`. Times leave out
  reading the input. Outputs of up to 10 qubits are checked against the
  input's unitary (up to a global phase); all pass.

`check.py` runs in CI on every pull request and push
(`.github/workflows/benchmarks.yml`): it compiles every circuit with QCC,
checks the small ones, and fails if any uses more gates than recorded in
`results/compilation.csv`, with Qiskit 2.6.0, which built the recorded
circuits. When QCC improves, rerun `run.py` to record it.

`plugin.py` compares Qiskit with and without QCC as its optimization stage
(below).

The suite follows the kinds of circuits of
[ucc-bench](https://github.com/unitaryfoundation/ucc-bench), but its code and
circuits are written here: ucc-bench is AGPL-licensed, this repository MIT.

## Results

Two-qubit gates / all gates (Intel i9-14900K, Qiskit 2.6.0, pytket 2.18.5,
UCC 0.4.12):

| circuit | qubits | QCC | Qiskit | pytket | UCC |
|---|---|---|---|---|---|
| QFT | 8 | 56 / 108 | 56 / 129 | 56 / 129 | 68 / 357 |
|  | 16 | 238 / 405 | 238 / 510 | 240 / 513 | 262 / 1285 |
|  | 32 | 686 / 1125 | 686 / 1454 | 962 / 1994 | 734 / 3593 |
| QAOA | 8 | 48 / 159 | 48 / 174 | 48 / 174 | 48 / 306 |
|  | 16 | 96 / 317 | 96 / 346 | 96 / 342 | 96 / 618 |
|  | 32 | 192 / 634 | 192 / 684 | 192 / 686 | 192 / 1239 |
| Quantum Volume | 8 | 72 / 472 | 72 / 472 | 72 / 731 | 72 / 753 |
|  | 16 | 369 / 2294 | 369 / 2294 | 369 / 3961 | 369 / 3708 |
|  | 32 | 1479 / 9034 | 1479 / 9034 | 1479 / 15749 | 1479 / 14822 |
| GHZ | 8 | 7 / 10 | 7 / 10 | 7 / 10 | 7 / 21 |
|  | 16 | 15 / 18 | 15 / 18 | 15 / 18 | 15 / 29 |
|  | 32 | 31 / 34 | 31 / 34 | 31 / 34 | 31 / 45 |
| Heisenberg | 8 | 84 / 295 | 84 / 356 | 84 / 434 | 84 / 674 |
|  | 16 | 180 / 627 | 180 / 760 | 180 / 914 | 180 / 1450 |
|  | 32 | 372 / 1291 | 372 / 1568 | 372 / 1874 | 372 / 3002 |
| Adder | 8 | 40 / 105 | 49 / 117 | 46 / 129 | 49 / 424 |
|  | 16 | 92 / 245 | 113 / 273 | 106 / 301 | 113 / 988 |
|  | 32 | 196 / 525 | 241 / 585 | 226 / 645 | 241 / 2116 |
| Clifford+T | 8 | 73 / 298 | 77 / 308 | 73 / 280 | 77 / 871 |
|  | 16 | 154 / 524 | 155 / 550 | 148 / 495 | 155 / 1743 |
|  | 32 | 313 / 1079 | 316 / 1104 | 310 / 1039 | 316 / 3450 |

Geometric means over the 21 circuits, relative to Qiskit, and median time:

| compiler | two-qubit gates | all gates | median time |
|---|---|---|---|
| QCC | 0.968× | 0.913× | 7 ms |
| Qiskit | 1.000× | 1.000× | 7 ms |
| pytket | 1.002× | 1.123× | 1161 ms |
| UCC | 1.017× | 2.210× | 21 ms |

- QCC uses the fewest gates of the four compilers: 3% fewer two-qubit gates
  than Qiskit (0.968×) and 9% fewer gates in all, and no more of either than
  Qiskit on any circuit. Phase folding (`core/src/phase_folding.rs`) does
  most of it: on the adder it removes 19% of the two-qubit gates (196 against
  Qiskit's 241 and pytket's 226 at 32 qubits), and on the QFT, QAOA and
  Heisenberg chain QCC leaves 7–23% fewer gates in all than Qiskit.
- Qiskit, pytket and QCC all remove the QFT's final SWAPs by relabelling the
  output qubits (3 CX per SWAP: 12, 24 and 48); `run.py` asks QCC to with
  `relabel_outputs=True` and checks its outputs with the permutation it
  reports. By default QCC only does so when the moved qubits are measured
  last, as an OpenQASM file cannot say that its qubits were relabelled.
- Where it may relabel, QCC's block resynthesis is also permutation-aware:
  a block U that is a SWAP away from needing fewer CX (a CX each way is a
  SWAP and one CX) is synthesized as SWAP·U, and the SWAP elided. On
  Clifford+T that saves 1–5% of the two-qubit gates, for a few more
  single-qubit gates; pytket, which does the same with its implicit wire
  swaps and more Clifford rewrites, still uses up to 4% fewer there.
- On Quantum Volume, QCC and Qiskit leave the same gates. Synthesizing each
  two-qubit block for its place in the circuit (`Objective::InCircuit` in
  `core/src/two_qubit.rs`: fewest gates between its first and last two-qubit
  gates, as the outer ones merge with their neighbours') saves 1% there over
  synthesizing it alone.
- QCC takes about as long as Qiskit on most circuits (a median of 7 ms for
  both), but up to 6× longer on Quantum Volume, where most of the work is
  resynthesizing generic two-qubit blocks.
- Like Qiskit's default, QCC drops two-qubit interactions whose removal
  costs less than 10⁻⁹ in average gate infidelity (the QFT's tiniest
  controlled phases); see `core/src/two_qubit.rs`.

## Inside Qiskit

With `qevedo-compiler` installed, Qiskit can use QCC as its optimization
stage: `transpile(circuit, backend, optimization_method="qcc")` keeps
Qiskit's layout, routing and translation and replaces its optimizations with
QCC's (`python/qevedo/compiler/qiskit_plugin.py`). `plugin.py` transpiles the
8- and 16-qubit circuits at optimization level 3 onto a 19-qubit heavy-hex
chip, with and without it, and checks every output's state
(`results/plugin.csv`). Geometric means of Qiskit + QCC over Qiskit alone:

| gates | two-qubit gates | all gates |
|---|---|---|
| {rz, sx, x, cx} | 0.973× | 0.867× |
| {rz, sx, x, cz} | 1.000× | 0.981× |

- With CX, QCC removes 13% of the gates and 3% of the two-qubit gates (the
  16-qubit QFT: 1324 → 518 gates; the adder: 180 → 159 two-qubit gates);
  Quantum Volume gets up to 1% more single-qubit gates.
- With CZ, QCC saves 2% of the gates but gives the Heisenberg chain 20–22%
  more single-qubit gates. Its two-qubit blocks there have a symmetric
  interaction (XX + YY + ZZ), whose decomposition leaves a choice of outer
  single-qubit gates that QCC does not yet make well; Qiskit picks simple
  ones. That is the next thing to fix.
