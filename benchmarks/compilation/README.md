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
| Heisenberg | 8 | 84 / 273 | 84 / 356 | 84 / 434 | 84 / 674 |
|  | 16 | 180 / 573 | 180 / 760 | 180 / 914 | 180 / 1450 |
|  | 32 | 372 / 1173 | 372 / 1568 | 372 / 1874 | 372 / 3002 |
| Adder | 8 | 40 / 105 | 49 / 117 | 46 / 129 | 49 / 424 |
|  | 16 | 92 / 245 | 113 / 273 | 106 / 301 | 113 / 988 |
|  | 32 | 196 / 525 | 241 / 585 | 226 / 645 | 241 / 2116 |
| Clifford+T | 8 | 73 / 293 | 77 / 308 | 73 / 280 | 77 / 871 |
|  | 16 | 154 / 525 | 155 / 550 | 148 / 495 | 155 / 1743 |
|  | 32 | 313 / 1074 | 316 / 1104 | 310 / 1039 | 316 / 3450 |
| Hamiltonian | 8 | 108 / 451 | 191 / 474 | 178 / 403 | 191 / 1986 |
|  | 16 | 325 / 1285 | 393 / 972 | 390 / 887 | 393 / 4186 |
|  | 32 | 740 / 1568 | 740 / 1789 | 737 / 1681 | 740 / 7739 |

Geometric means over the 21 circuits, relative to Qiskit, and median time:

| compiler | two-qubit gates | all gates | median time |
|---|---|---|---|
| QCC | 0.942× | 0.916× | 4 ms |
| Qiskit | 1.000× | 1.000× | 6 ms |
| pytket | 0.998× | 1.093× | 978 ms |
| UCC | 1.015× | 2.400× | 14 ms |

- QCC uses the fewest gates of the four compilers: 6% fewer two-qubit gates
  than Qiskit (0.942×) and 8% fewer gates in all, and no more two-qubit
  gates than Qiskit on any circuit. Phase folding
  (`core/src/phase_folding.rs`) does most of it on the arithmetic: on the
  adder it removes 19% of the two-qubit gates (196 against Qiskit's 241 and
  pytket's 226 at 32 qubits), and on the QFT, QAOA and Heisenberg chain QCC
  leaves 7–25% fewer gates in all than Qiskit.
- On one Trotter step of a random Hamiltonian (`hamiltonian` in
  `circuits.py`: Pauli strings of weight 2–6, with Qiskit's CX ladders), the
  Pauli network (`core/src/pauli_network.rs`) builds the circuit again with
  shared CNOTs: 108 against Qiskit's 191 at 8 qubits, and 325 against 393 at
  16, where it leaves more single-qubit gates (1285 against 972) — the one
  circuit where QCC uses more gates in all. At 32 qubits it does not beat
  the direct compilation, which QCC then keeps.
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
- On the Heisenberg chain, moving single-qubit gates through its symmetric
  XX + YY + ZZ blocks (`core/src/symmetric_blocks.rs`) saves another 8–9%
  of the gates: 1173 at 32 qubits, against Qiskit's 1568.
- QCC is faster than Qiskit on most circuits (a median of 4 ms, against
  6), but 2–3× slower on Quantum Volume, where most of the work is
  resynthesizing generic two-qubit blocks (on all cores), and up to 10×
  on the Hamiltonian, where it also tries the Pauli network.
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
| {rz, sx, x, cx} | 0.976× | 0.849× |
| {rz, sx, x, cz} | 0.978× | 0.895× |

- With CX, QCC removes 15% of the gates and 2% of the two-qubit gates (the
  16-qubit QFT: 1324 → 518 gates; the adder: 180 → 159 two-qubit gates);
  Quantum Volume gets up to 1% more single-qubit gates.
- With CZ, QCC removes 10% of the gates, and 7–13% on the Heisenberg chain
  (16 qubits: 992 → 866). Its XX + YY + ZZ blocks have a symmetric Weyl
  point, where the outer single-qubit gates of the KAK decomposition are
  not unique; QCC now picks them for each block and its neighbours
  (`symmetric_variants` in `core/src/two_qubit.rs` and
  `core/src/symmetric_blocks.rs`). Before, it gave this circuit a fifth more
  gates than Qiskit.

## On a chip

`device.py` compiles the 8- and 16-qubit circuits onto the 19-qubit
heavy-hex chip with QCC alone (`compile_circuit` with the chip's couplings:
it optimizes, lays out, routes and optimizes again; see
`docs/ARCHITECTURE.md`), with Qiskit at level 3, and with Qiskit using QCC as
its optimization stage, and checks every output's state
(`results/device.csv`). Geometric means over Qiskit's level 3:

| gates | compiler | two-qubit gates | all gates | median time |
|---|---|---|---|---|
| {rz, sx, x, cx} | QCC | 0.943× | 0.835× | 50 ms |
|  | Qiskit + QCC | 0.976× | 0.849× | 45 ms |
| {rz, sx, x, cz} | QCC | 0.937× | 0.881× | 69 ms |
|  | Qiskit + QCC | 0.978× | 0.895× | 96 ms |

Two-qubit gates / all gates with CX:

| circuit | qubits | QCC | Qiskit | Qiskit + QCC |
|---|---|---|---|---|
| QFT | 8 | 84 / 136 | 102 / 188 | 98 / 150 |
|  | 16 | 357 / 524 | 356 / 1324 | 350 / 518 |
| QAOA | 8 | 94 / 207 | 97 / 232 | 97 / 209 |
|  | 16 | 188 / 411 | 187 / 445 | 187 / 408 |
| Quantum Volume | 8 | 111 / 515 | 108 / 575 | 108 / 560 |
|  | 16 | 822 / 2787 | 828 / 3182 | 828 / 3203 |
| GHZ | 8 | 7 / 10 | 7 / 10 | 7 / 10 |
|  | 16 | 15 / 18 | 15 / 18 | 15 / 18 |
| Heisenberg | 8 | 84 / 271 | 84 / 356 | 84 / 300 |
|  | 16 | 180 / 571 | 180 / 760 | 180 / 631 |
| Adder | 8 | 57 / 151 | 72 / 184 | 59 / 161 |
|  | 16 | 148 / 405 | 180 / 452 | 159 / 409 |
| Clifford+T | 8 | 236 / 521 | 249 / 546 | 249 / 484 |
|  | 16 | 751 / 1279 | 774 / 1394 | 772 / 1211 |
| Hamiltonian | 8 | 352 / 796 | 425 / 867 | 425 / 693 |
|  | 16 | 1403 / 2113 | 1463 / 2348 | 1463 / 1995 |

- QCC's own layout and routing do better than Qiskit's on the adder (148
  CNOTs against 180) and Clifford+T, and as well on the rest, give or take a
  few gates. GHZ and the Heisenberg chain fit the chip, and both find the
  layout that needs no SWAP.
- Its SABRE counts a SWAP right after a block on the same pair for the
  CNOTs it costs once merged into the block, as NASSC does (Liu, Li and
  Zhou, HPCA 2022). That, and 32 layout trials instead of 8, took 8% of the
  two-qubit gates off: on the 16-qubit QFT, 479 CNOTs became 357, against
  Qiskit's 356.
- QCC uses a few more two-qubit gates than Qiskit on three circuits:
  Quantum Volume at 8 qubits (111 against 108, with CX or CZ), QAOA at 16
  with CX (188 against 187) and the QFT at 16 (357 against 356).
- On the Hamiltonian, a Pauli network can need more routing than it saves,
  so QCC routes its logical circuit with one and without, and keeps the
  cheaper: 1403 CNOTs against Qiskit's 1463 at 16 qubits.
