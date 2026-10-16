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
- `run.py` compiles them with QCC's default pipeline, Qiskit's `transpile`
  at optimization level 3, pytket's `FullPeepholeOptimise` plus a rebase,
  and UCC's `compile`, and writes `results/compilation.csv`. Times leave out
  reading the input. Outputs of up to 10 qubits are checked against the
  input's unitary (up to a global phase); all pass.

The suite follows the kinds of circuits of
[ucc-bench](https://github.com/unitaryfoundation/ucc-bench), but its code and
circuits are written here: ucc-bench is AGPL-licensed, this repository MIT.

## Results

Two-qubit gates / all gates (Intel i9-14900K, Qiskit 2.6.0, pytket 2.18.5,
UCC 0.4.12):

| circuit | qubits | QCC | Qiskit | pytket | UCC |
|---|---|---|---|---|---|
| QFT | 8 | 68 / 141 | 56 / 129 | 56 / 129 | 68 / 357 |
| | 16 | 262 / 534 | 238 / 510 | 240 / 513 | 262 / 1285 |
| | 32 | 734 / 1502 | 686 / 1454 | 962 / 1994 | 734 / 3593 |
| QAOA | 8 | 48 / 174 | 48 / 174 | 48 / 174 | 48 / 306 |
| | 16 | 96 / 346 | 96 / 346 | 96 / 342 | 96 / 618 |
| | 32 | 192 / 684 | 192 / 684 | 192 / 686 | 192 / 1239 |
| Quantum Volume | 8 | 72 / 528 | 72 / 472 | 72 / 731 | 72 / 753 |
| | 16 | 369 / 2575 | 369 / 2294 | 369 / 3961 | 369 / 3708 |
| | 32 | 1479 / 10527 | 1479 / 9034 | 1479 / 15749 | 1479 / 14822 |
| GHZ | 8 | 7 / 10 | 7 / 10 | 7 / 10 | 7 / 21 |
| | 16 | 15 / 18 | 15 / 18 | 15 / 18 | 15 / 29 |
| | 32 | 31 / 34 | 31 / 34 | 31 / 34 | 31 / 45 |
| Heisenberg | 8 | 84 / 377 | 84 / 356 | 84 / 434 | 84 / 674 |
| | 16 | 180 / 761 | 180 / 760 | 180 / 914 | 180 / 1450 |
| | 32 | 372 / 1529 | 372 / 1568 | 372 / 1874 | 372 / 3002 |
| Adder | 8 | 49 / 117 | 49 / 117 | 46 / 129 | 49 / 424 |
| | 16 | 113 / 273 | 113 / 273 | 106 / 301 | 113 / 988 |
| | 32 | 241 / 585 | 241 / 585 | 226 / 645 | 241 / 2116 |
| Clifford+T | 8 | 77 / 290 | 77 / 308 | 73 / 280 | 77 / 871 |
| | 16 | 155 / 522 | 155 / 550 | 148 / 495 | 155 / 1743 |
| | 32 | 316 / 1070 | 316 / 1104 | 310 / 1039 | 316 / 3450 |

Geometric means over the 21 circuits, relative to Qiskit, and median time:

| compiler | two-qubit gates | all gates | median time |
|---|---|---|---|
| QCC | 1.017× | 1.021× | 2 ms |
| Qiskit | 1.000× | 1.000× | 9 ms |
| pytket | 1.002× | 1.123× | 875 ms |
| UCC | 1.017× | 2.210× | 16 ms |

- QCC uses as many two-qubit gates as Qiskit on every circuit but the QFT,
  where Qiskit and pytket remove the final SWAPs by relabelling the output
  qubits (3 CX per SWAP: 12, 24 and 48). QCC keeps them, since its OpenQASM
  output cannot say that the qubits were relabelled.
- pytket finds fewer two-qubit gates on the adder and Clifford+T circuits
  (up to 6%), at a few hundred times the time.
- QCC needs fewer gates in all than Qiskit on Clifford+T and on the
  32-qubit Heisenberg chain, but 12–17% more on Quantum Volume, all of them
  single-qubit gates.
- QCC is faster than Qiskit on 16 of the 21 circuits, but slower on Quantum
  Volume (0.32 s against 0.07 s at 32 qubits) and the Heisenberg chain, where
  most of the work is resynthesizing generic two-qubit blocks.
- Like Qiskit's default, QCC drops two-qubit interactions whose removal
  costs less than 10⁻⁹ in average gate infidelity (the QFT's tiniest
  controlled phases); see `core/src/two_qubit.rs`.
