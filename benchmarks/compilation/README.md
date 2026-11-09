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

Two-qubit gates / all gates (a 4-vCPU Intel Xeon VM at 2.8 GHz, Qiskit 2.6.0,
pytket 2.18.5, UCC 0.4.12):

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
| Clifford+T | 8 | 73 / 278 | 77 / 308 | 73 / 280 | 77 / 871 |
|  | 16 | 150 / 505 | 155 / 550 | 148 / 495 | 155 / 1743 |
|  | 32 | 310 / 1034 | 316 / 1104 | 310 / 1039 | 316 / 3450 |
| Hamiltonian | 8 | 101 / 286 | 191 / 474 | 178 / 403 | 191 / 1986 |
|  | 16 | 313 / 790 | 393 / 972 | 390 / 887 | 393 / 4186 |
|  | 32 | 687 / 1511 | 740 / 1789 | 737 / 1681 | 740 / 7739 |

Geometric means over the 21 circuits, relative to Qiskit, and median time:

| compiler | two-qubit gates | all gates | median time |
|---|---|---|---|
| QCC | 0.933× | 0.875× | 7 ms |
| Qiskit | 1.000× | 1.000× | 11 ms |
| pytket | 0.998× | 1.093× | 1803 ms |
| UCC | 1.015× | 2.400× | 33 ms |

- QCC uses the fewest gates of the four compilers: 7% fewer two-qubit gates
  than Qiskit (0.933×) and 12% fewer gates in all, and no more two-qubit
  gates than Qiskit on any circuit. Phase folding
  (`core/src/phase_folding.rs`) does most of it on the arithmetic: on the
  adder it removes 19% of the two-qubit gates (196 against Qiskit's 241 and
  pytket's 226 at 32 qubits), and on the QFT, QAOA and Heisenberg chain QCC
  leaves 7–25% fewer gates in all than Qiskit.
- On one Trotter step of a random Hamiltonian (`hamiltonian` in
  `circuits.py`: Pauli strings of weight 2–6, with Qiskit's CX ladders), the
  Pauli network (`core/src/pauli_network.rs`) builds the circuit again with
  shared CNOTs: 101 against Qiskit's 191 at 8 qubits, 313 against 393 at
  16 and 687 against 740 at 32, with fewer gates in all too (286 against
  474, 790 against 972 and 1511 against 1789; pytket leaves 887 at 16 and
  737 / 1681 at 32). Single-qubit gates used to be its weakness (1285
  in all at 16 qubits, more than Qiskit): the local Cliffords that set up
  each CX are no longer undone after it, each CX is placed whichever way
  round needs fewer of them, basis changes use √X rather than H, and eight
  syntheses with differently seeded tie-breaking (as pytket's
  GreedyPauliSimp offers with `trials`) run side by side (twelve since
  0.4.7, below). On 15 other random
  Hamiltonians of 6 to 24 qubits this left a geometric mean of 0.934× the
  CNOTs and 0.779× the gates of 0.4.0, and more of neither on any.
- QCC keeps a network only when it needs no more two-qubit gates and no
  more gates in all than the direct compilation, and judges each trial on
  its own (0.4.7): at 32 qubits the trial with the fewest CNOTs (663) needs
  1581 gates in all, more than the direct compilation's 1555, and before
  0.4.7 that trial was the only one judged, so the direct compilation was
  kept (740 / 1555). The greedy synthesis also charges each local Clifford
  a candidate CX needs; with the usual charge, a hundredth of a weight, it
  only breaks ties. Four more trials (seeds 0–3) charge a tenth of a
  weight, giving up a few CNOTs for fewer single-qubit gates: at 32 qubits
  one of them gives 687 / 1511, which is kept; at 8 qubits 101 / 286
  (against 107 / 293); at 16 the usual trials stay best (the larger charge
  finds 335 CNOTs there). The eight usual trials still run, so a circuit
  can only get cheaper. On 15 random Hamiltonians of 6–24 qubits not used
  for tuning (seeds other than the benchmark's), this left 0.981× the
  CNOTs and 0.980× the gates of 0.4.6, and more of neither on any; it costs
  1.3–1.5× the time on the Hamiltonians (32 qubits: 0.33 → 0.40–0.45 s,
  medians of five, side by side), and nothing on the other circuits, where
  no network is tried.
- Qiskit, pytket and QCC all remove the QFT's final SWAPs by relabelling the
  output qubits (3 CX per SWAP: 12, 24 and 48); `run.py` asks QCC to with
  `relabel_outputs=True` and checks its outputs with the permutation it
  reports. By default QCC only does so when the moved qubits are measured
  last, as an OpenQASM file cannot say that its qubits were relabelled.
- Where it may relabel, QCC's block resynthesis is also permutation-aware:
  a block U that is a SWAP away from needing fewer CX (a CX each way is a
  SWAP and one CX) is synthesized as SWAP·U, and the SWAP elided.
- Blocks on the same pair that a gate on one of their qubits and a third
  keeps apart are gathered when that gate commutes out of the way
  (`core/src/commuting_blocks.rs`, as tket's Clifford simplification does
  for CXs), and single-qubit gates are traded through each lone CX, a Z
  rotation on its control and an X rotation on its target, to make the
  runs around it cheapest (`core/src/symmetric_blocks.rs`). On Clifford+T
  that took the 16-qubit circuit from 154 / 525 to 150 / 505 and the
  32-qubit one from 313 / 1074 to 310 / 1034: as many CNOTs as pytket at 32
  qubits and fewer gates in all (1034 against 1039); at 16 pytket still
  uses 2 fewer CNOTs (148). The Hamiltonian lost 2–13 gates (792 → 790 at
  16 qubits, 1568 → 1555 at 32).
- On the adder, gathering would leave 37, 85 and 181 CNOTs instead of 40,
  92 and 196, but 108, 252 and 540 gates in all instead of 105, 245 and
  525: the SWAPs elided from the gathered blocks leave generic single-qubit
  gates behind (pytket's Clifford simplification, run on QCC's output,
  makes the same trade: 37 / 115 at 8 qubits). QCC keeps a gathered circuit
  only when it needs no more two-qubit gates and no more gates in all, so
  the adder is unchanged here (on the chip with CZ, where the trade pays,
  the 8-qubit adder went from 57 / 235 to 54 / 230).
- On Quantum Volume, QCC and Qiskit leave the same gates. Synthesizing each
  two-qubit block for its place in the circuit (`Objective::InCircuit` in
  `core/src/two_qubit.rs`: fewest gates between its first and last two-qubit
  gates, as the outer ones merge with their neighbours') saves 1% there over
  synthesizing it alone.
- On the Heisenberg chain, moving single-qubit gates through its symmetric
  XX + YY + ZZ blocks (`core/src/symmetric_blocks.rs`) saves another 8–9%
  of the gates: 1173 at 32 qubits, against Qiskit's 1568.
- QCC is faster than Qiskit on most circuits (18 of 24; a median of 7 ms,
  against 11), but slower on the Hamiltonian (3.8–27× in this run), where
  it also tries the Pauli networks (twelve of them, side by side), on
  Quantum Volume at 16 and 32 qubits (1.8× and 5.0×), where most of the
  work is resynthesizing generic two-qubit blocks (on all cores, of which
  this machine has four), and by a fraction of a millisecond on Clifford+T
  at 32 (17.4 against 16.7 ms). Times are single compilations and vary
  between runs on this VM: in an earlier run the 32-qubit Quantum Volume
  circuit's 316 ms was a median of 191 ms over nine compilations, against
  Qiskit's 69 ms.
- 0.4.3 made every circuit faster without changing a gate (the outputs are
  identical): a median of 11 compilations of each circuit, 0.4.2 and 0.4.3
  side by side on this machine, gave 0.73× the time in geometric mean, and
  0.56–0.92× on each; Quantum Volume at 32 qubits went from 230 to 155 ms,
  the Rust core's share from about 165 to 95 ms. Counting single-qubit gates
  (what two-qubit synthesis does most, for every way of moving a rotation
  across a gate) now reads the Euler angles' tests off the unitary's entries
  instead of computing the angles with `atan2` (`Shape` in
  `core/src/one_qubit.rs`); phase folding keeps parities as short lists of
  variables rather than bitsets over thousands of them; rotation merging
  updates its Clifford's images on the qubits a gate touches only; and block
  resynthesis works out each block's unitary and CX counts once instead of
  three times. Of the rest, about a third of the time at 32 qubits is the
  Python wrapper around the core (converting and copying the instructions).
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
| {rz, sx, x, cx} | 0.970× | 0.844× |
| {rz, sx, x, cz} | 0.936× | 0.837× |

- With CX, QCC removes 16% of the gates and 3% of the two-qubit gates (the
  16-qubit QFT: 1324 → 516 gates; the adder: 180 → 159 two-qubit gates);
  Quantum Volume gets up to 1% more single-qubit gates.
- With CZ, QCC removes 16% of the gates and 6% of the two-qubit gates
  (the 16-qubit QFT: 380 → 348 CZs, 2022 → 1244 gates; the 16-qubit
  Hamiltonian: 1477 → 1312 CZs), and 15–16% of the gates on the Heisenberg
  chain (16 qubits: 992 → 844). Its XX + YY + ZZ blocks have a symmetric Weyl
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
(`results/device.csv`; last run for 0.4.8, on a 4-vCPU Xeon VM). Geometric means over
Qiskit's level 3:

| gates | compiler | two-qubit gates | all gates | median time |
|---|---|---|---|---|
| {rz, sx, x, cx} | QCC | 0.928× | 0.818× | 32 ms |
|  | Qiskit + QCC | 0.970× | 0.844× | 43 ms |
| {rz, sx, x, cz} | QCC | 0.899× | 0.834× | 111 ms |
|  | Qiskit + QCC | 0.936× | 0.837× | 142 ms |

Qiskit's own median in the same run: 34 ms with CX, 36 ms with CZ.

Two-qubit gates / all gates with CX:

| circuit | qubits | QCC | Qiskit | Qiskit + QCC |
|---|---|---|---|---|
| QFT | 8 | 80 / 132 | 102 / 188 | 92 / 144 |
|  | 16 | 353 / 520 | 356 / 1324 | 348 / 516 |
| QAOA | 8 | 94 / 207 | 97 / 232 | 95 / 207 |
|  | 16 | 179 / 403 | 187 / 445 | 183 / 405 |
| Quantum Volume | 8 | 111 / 515 | 108 / 575 | 108 / 560 |
|  | 16 | 783 / 2781 | 828 / 3182 | 828 / 3203 |
| GHZ | 8 | 7 / 10 | 7 / 10 | 7 / 10 |
|  | 16 | 15 / 18 | 15 / 18 | 15 / 18 |
| Heisenberg | 8 | 84 / 271 | 84 / 356 | 84 / 300 |
|  | 16 | 180 / 571 | 180 / 760 | 180 / 631 |
| Adder | 8 | 57 / 151 | 72 / 184 | 59 / 161 |
|  | 16 | 148 / 405 | 180 / 452 | 159 / 409 |
| Clifford+T | 8 | 236 / 512 | 249 / 546 | 249 / 479 |
|  | 16 | 733 / 1209 | 774 / 1394 | 772 / 1205 |
| Hamiltonian | 8 | 345 / 648 | 425 / 867 | 425 / 686 |
|  | 16 | 1336 / 2100 | 1463 / 2348 | 1463 / 1985 |

- QCC's own layout and routing do better than Qiskit's on the adder (148
  CNOTs against 180) and Clifford+T, and as well on the rest, give or take a
  few gates. GHZ and the Heisenberg chain fit the chip, and both find the
  layout that needs no SWAP.
- Its SABRE counts a SWAP right after a block on the same pair for the
  CNOTs it costs once merged into the block, as NASSC does (Liu, Li and
  Zhou, HPCA 2022). That, and 32 layout trials instead of 8, took 8% of the
  two-qubit gates off: on the 16-qubit QFT, 479 CNOTs became 357.
- Gathering blocks through commuting gates (`core/src/commuting_blocks.rs`)
  matters more after routing, where SWAPs sit next to blocks on other
  pairs: with CX it took the 16-qubit QFT from 357 to 353 CNOTs (Qiskit
  356), QAOA from 188 to 180 (Qiskit 187), the 16-qubit Hamiltonian from
  1403 to 1368, and with CZ the 16-qubit Clifford+T circuit from 765 to 722
  CZs (2.4–7% of the two-qubit gates on each). Inside Qiskit it took the
  QFT at 16 qubits from 380 to 348 CZs. Four rows got worse in one count or
  both, as the circuits that reach routing differ: with CX, Clifford+T at
  16 qubits (752 / 1251 against 751 / 1279); with CZ, the adder at 16 (150 /
  616 against 148 / 586), the Hamiltonian at 8 (336 / 1378 against 337 /
  1339) and the QFT at 16 (353 / 1422 against 357 / 1408). Since the
  repeated rounds above, the Hamiltonian at 8 is better on both (331 /
  1333, and 310 / 1119 since 0.4.7); the adder (149 / 615) and the QFT (353 / 1421) still have more
  gates in all than before gathering.
- For CZ the pipeline is not a fixed point, as it nearly is for CX: it
  reads the circuit through CX (each CZ as H·CX·H) so that phase folding
  sees parities, then translates it back, and each translation leaves the
  single-qubit gates where the next reading finds parities and blocks the
  last did not. QCC therefore runs the pipeline again on its own result,
  on the physical qubits and inside Qiskit, while that gets cheaper, up to
  eight times (`until_fixed_point` in `core/src/passes.rs`), as Qiskit's
  level 3 repeats its optimization loop. With CZ, every circuit got cheaper
  or stayed the same, none got worse: QCC alone 0.917× → 0.908× the
  two-qubit gates of Qiskit and 0.868× → 0.848× the gates; the 16-qubit
  Hamiltonian 1379 / 5372 → 1292 / 4572 (Qiskit 1477 / 5357), Clifford+T
  at 16 qubits 722 / 2880 → 708 / 2790, so QCC no longer uses more gates
  in all than Qiskit there (2844). Inside Qiskit: 0.947× / 0.861× →
  0.936× / 0.837×; Clifford+T at 16 qubits 756 / 2579 → 741 / 2455. It
  costs time with CZ only: on this machine, against the previous version
  side by side (median of five), QCC alone took 1.65× as long in geometric
  mean over the CZ circuits (the 16-qubit Hamiltonian 1.7 s → 4.4 s, where
  Qiskit takes 0.13 s), and the pipeline inside Qiskit 2.5× (tens of
  milliseconds). With CX nothing changes. 0.4.5 won back most of that
  (below).
- 0.4.5 runs work that reads only the same input side by side, on all
  cores (`core/src/parallel.rs`): SABRE's 32 layout trials, the direct
  and through-CX compilations for CZ (in every round above), the Pauli
  networks beside them, and the routing of the logical circuit with and
  without a network, which before compiled the direct candidates twice.
  Results come back in the order a loop would give them and the first of
  the cheapest is kept, so the outputs are byte-identical: all 112
  compilations of the benchmark circuits (all-to-all, on the chip and
  inside Qiskit, with CX and CZ) give the same QASM as 0.4.4. Time on this
  4-core machine, 0.4.4 and 0.4.5 side by side (best of two medians of
  five): on the chip 0.51× with CX and 0.62× with CZ in geometric mean
  (the 16-qubit Hamiltonian 1.38 s → 0.57 s with CX, 5.58 s → 2.21 s with
  CZ; Quantum Volume 16 0.46 s → 0.18 s with CX); inside Qiskit with CZ
  0.82×, all-to-all with CZ 0.85×. Where nothing new runs side by side
  (CX all-to-all and inside Qiskit) the geometric means were 1.04× and
  1.00× (0.98× in an earlier run): noise, between 0.8× and 1.4× on
  circuits of a few milliseconds. Circuits under 64 instructions are
  compiled on one thread, as a thread costs as much as their compilation.
  QCC is still slower than Qiskit on the CZ chip's largest circuits (the
  16-qubit Hamiltonian 2.2 s against 0.15 s), where the rounds remain.
- 0.4.6 makes each round cheaper, again with byte-identical outputs (the
  same 112 compilations). All the pipelines of one compilation share
  their two-qubit syntheses (`SynthesisCache` in `core/src/lower.rs`): on
  the 16-qubit Hamiltonian on the chip with CZ, 17 374 of the 26 156
  synthesis requests (66%) had been synthesized before by another
  pipeline or round of the same compilation, 4 341 of 6 251 on
  Clifford+T. And gathering blocks (`core/src/commuting_blocks.rs`) no
  longer copies the whole region between two blocks for every cut it
  tries, nor works out the same range's CX count again for each cut:
  17–145 ms → 4–45 ms per pipeline there. Against 0.4.5 side by side
  (best of two medians of five): on the chip with CZ 0.83× in geometric
  mean (the 16-qubit Hamiltonian 2.14 s → 1.60 s, Clifford+T 16 0.72 s →
  0.46 s), inside Qiskit with CZ 0.82× (the Hamiltonian 1.38 s →
  0.97 s), all-to-all with CZ 0.94×; with CX 0.96×, 1.00× and 1.02× (noise,
  0.8–1.4× on circuits of a few milliseconds). The 16-qubit Hamiltonian
  with CZ still takes 1.6 s against Qiskit's 0.15 s.
- QCC now uses more two-qubit gates than Qiskit on one circuit: Quantum
  Volume at 8 qubits (111 against 108, with CX or CZ). That row is SABRE's
  luck of the draw more than a gap: Qiskit's level 3 leaves 108 or 111
  there depending on its seed (seeds 1–10: four 108s, six 111s), and QCC's
  estimate of what its SWAPs cost is exact on it (72 CNOTs before routing,
  39 for the SWAPs), but no wider search tried (128 layouts, or 32
  routings each) found a routing for less.
- 0.4.8 routes each layout 16 times with different tie-breaking instead of
  4, and also optimizes the routing that the first 4 would have given,
  keeping the cheaper: the SWAPs' CNOTs are only an estimate of the
  optimized circuit, and on 40 held-out circuits (Quantum Volume, Clifford+T,
  Hamiltonians and QAOA of 10–16 qubits, other seeds than these) the wider
  search alone left 3 with more two-qubit gates. With both, 19 of the 40
  got fewer two-qubit gates and none more (19 213 → 18 936 in all, and
  66 134 → 64 817 gates); three trade a few single-qubit gates for fewer
  two-qubit ones, as QCC's cost ordering does everywhere. On the benchmark,
  seven rows got cheaper and none dearer: Quantum Volume at 16 qubits 822
  → 783 two-qubit gates with CX or CZ (Qiskit 828; Qiskit's own seeds give
  798–840 there), Clifford+T 752 / 1251 → 733 / 1209 with CX and 708 /
  2790 → 692 / 2677 with CZ, the Hamiltonian 1368 / 2106 → 1336 / 2100
  with CX and 1292 / 4572 → 1256 / 4454 with CZ, QAOA 180 / 409 → 179 /
  403 with CX. SABRE itself got cheaper in the same change, making the
  same choices (`core/src/routing/sabre.rs`; 2.1–2.9× fewer instructions
  for a search on four 16-qubit circuits), so the wider search costs
  little: against 0.4.7 side by side (best of two medians of five), 0.86×
  the time with CX and 0.97× with CZ in geometric mean. Where both
  routings go through the CZ rounds it costs more: the 16-qubit
  Hamiltonian with CZ 1.70 s → 2.33 s (1.37×), Quantum Volume at 16 with
  CZ 1.22×.
- On the Hamiltonian, a Pauli network can need more routing than it saves,
  so QCC routes its logical circuit with one and without, and keeps the
  cheaper: 345 CNOTs against Qiskit's 425 at 8 qubits (648 gates in all,
  against 867), and 1368 against 1463 at 16. The trials with the larger
  local-Clifford charge (0.4.7, above) took the 8-qubit one from 349 / 713
  to 345 / 648 with CX, and from 331 / 1333 to 310 / 1119 CZs / gates with
  CZ (Qiskit 427 / 1720); no other circuit on the chip changed.
