# Improvement log

One entry per change made by the daily improvement routine: what changed,
why, what it measured, and what to try next. Numbers come from the scripts
in `benchmarks/`.

## 2026-11-01 — Faster compilation, the same gates

**What.** Four exact speed-ups in the Rust core, none changing a gate:

- `core/src/one_qubit.rs`: `count_1q`, the single-qubit gate count that
  two-qubit synthesis calls for every way of moving a rotation across a gate
  (some 650 000 times on the 32-qubit Quantum Volume circuit), no longer
  computes the ZYZ angles. A count only asks whether certain angles are
  multiples of 2π, and those angles are the phases of products of the
  unitary's entries (`e^{i(φ+λ)} ∝ u₁₁² conj(det u)`, `e^{iλ} ∝ u₁₁ conj(u₁₀)`,
  …), which `Shape` tests with products and square roots; θ is compared
  through its tangent or cosine. Where `euler_zyz` drops the phase of
  exactly one tiny entry, φ and λ alone need the square root's sign, and the
  count falls back to the angles.
- `core/src/phase_folding.rs`: parities are sorted lists of variables, not
  bitsets over every variable created (thousands on Quantum Volume, where
  each rotation's key was tens of words to copy and hash).
- `core/src/pauli_network.rs`: rotation merging updates its Clifford's images
  by visiting only the qubits in a generator's support, and a single-qubit
  gate's rows once instead of twice.
- `core/src/passes.rs`: block resynthesis works out each block's unitary and
  whether SWAP·U is worth synthesizing once (`Plan`), instead of three times
  (for the requests and for each measure), and no longer copies blocks.

**Why.** Quantum Volume compile time was the oldest open problem in this
log (named under "Next" twice): 3.2× Qiskit's time at 32 qubits on this
machine (median of five: 250 ms against 77). A profile (callgrind) showed
`atan2`, `hypot` and the complex square root inside `euler_zyz` at a quarter
of all instructions, and rotation merging, phase folding's hashing and
repeated KAK decompositions at most of the serial rest. These are exact
changes whose effect can be checked bit for bit, which a new synthesis
method could not offer.

**Measured** (4-vCPU Xeon VM). Outputs: the 112 compilations of the
benchmark circuits (all-to-all with CX and CZ, the heavy-hex chip with CX and
CZ, inside Qiskit) give byte-identical QASM before and after; `run.py`,
`device.py` and `plugin.py`, complete runs, record the same gate counts.
Time, 0.4.2 and this version side by side, median of 11 compilations per
circuit: every circuit faster, 0.73× in geometric mean (0.56–0.92×); Quantum
Volume 9.3 / 45.0 / 229.5 ms → 7.5 / 29.6 / 154.8 ms at 8 / 16 / 32 qubits;
the Hamiltonian 38 / 132 / 336 → 31 / 121 / 281 ms. The Rust core alone on
Quantum Volume 32: about 165 → 95 ms. Qiskit takes 69 ms there (median of
nine), so QCC is still 2.8× slower on it, a third of that now in the Python
wrapper. The recorded single-run medians: 7 ms against Qiskit's 16.

**Learned.** With four cores the parallel block syntheses were not the
bottleneck: the serial passes were (rotation merging alone was a quarter of
the main thread). A count test at exactly the tolerance (θ − π/2 = 1e-10) can
go either way under rounding in both methods; the new test checks to either
side of it.

**Next.**
1. The Python wrapper: about a third of Quantum Volume 32's time is
   converting and copying instructions (`PassManager.run` copies the circuit,
   `NativePipeline` builds 9 000 dataclasses).
2. Win back the adder's gathered CNOTs without the single-qubit cost (see
   2026-10-30).
3. The Hamiltonian's time (8 Pauli-network trials, 2.4–12× Qiskit): cache
   the templates' KAK per Weyl point, or stop trials early once one is
   clearly worse.

## 2026-10-30 — Gathering blocks through commuting gates, and single-qubit trades through CXs

**What.**
- A new pass, `core/src/commuting_blocks.rs`, before block resynthesis:
  when two two-qubit gates on the same pair are kept apart only by
  two-qubit gates that touch one of their qubits (at most three), those
  gates are moved out of the way, to before the first or after the second,
  or further into the runs around them, wherever the cut leaves the fewest
  CX-like gates. A gate moves only if it commutes exactly with the product
  of the pair's gates it crosses (an 8×8 check, so single-qubit gates count
  together), carrying with it the gates on its third qubit that depend on
  it, and only when the gathered blocks need fewer CX-like gates (counting
  SWAP·U when relabelling is allowed).
- The pipeline finishes the circuit both with and without gathering
  whenever gathering moved a gate, and keeps the gathered one only when it
  needs no more two-qubit gates and no more gates in all (the Pauli
  network's rule).
- `core/src/symmetric_blocks.rs`: a block at the point (a, 0, 0), a lone CX
  for one, trades an X rotation in its canonical frame on each qubit
  separately (a Z rotation on a CX's control, an X rotation on its target),
  chosen to make the runs before and after cheapest. As it turns Z
  rotations into generic runs that phase folding cannot merge, it runs only
  in the last pipeline of a flow (`last` in `best_pipeline`): not before
  routing, nor in the CX stage of compiling for CZ.

**Why.** pytket's Clifford simplification, run on QCC's own output,
still removed CNOTs: on the adder 40 → 37, 92 → 85, 196 → 181, on
Clifford+T 154 → 148. Tracing the smallest window where it did showed
pairs of CXs on one pair kept apart by a CX that commutes with one of
them: QCC's blocks are maximal runs as they stand in the list, so its
resynthesis never saw them together. The technique is published (Fagan
and Duncan, arXiv:1901.10114; tket, Sivarajah et al., arXiv:2003.10611),
is exact, and could be verified end to end today; the other candidates
(Quantum Volume speed, spread over every pass; Pauli-network gates)
promised less.

**Measured** (4-vCPU Xeon VM; `run.py`, `device.py`, `plugin.py`, complete
runs; `check.py` passes, no circuit worse in either count). All-to-all:
Clifford+T 73 / 293 → 73 / 278, 154 / 525 → 150 / 505, 313 / 1074 →
310 / 1034 (pytket 310 / 1039 at 32 qubits, 148 / 495 at 16); Hamiltonian
296 → 293, 792 → 790, 1568 → 1555 gates. Geometric means over Qiskit:
0.938× two-qubit gates (was 0.940×), 0.877× all gates (was 0.882×). On the
heavy-hex chip with CX: 0.935× / 0.825× (was 0.942× / 0.830×); the QFT at
16 qubits 357 → 353 CNOTs (Qiskit 356), QAOA 188 → 180 (Qiskit 187), the
Hamiltonian at 16 1403 → 1368; with CZ 0.917× / 0.868× (was 0.935× /
0.878×). Inside Qiskit: 0.970× / 0.844× with CX (was 0.976× / 0.849×),
0.947× / 0.861× with CZ (was 0.978× / 0.895×). Worse in one count or
both on the chip: CX Clifford+T 16 (752 / 1251 against 751 / 1279); CZ
adder 16 (150 / 616 against 148 / 586), Hamiltonian 8 (336 / 1378 against
337 / 1339), QFT 16 (353 / 1422 against 357 / 1408); inside Qiskit with
CZ, Clifford+T 8 (245 / 885 against 251 / 847) and QAOA 8 (95 / 367
against 97 / 366). Time: Quantum Volume and Heisenberg unchanged; circuits
where gathering moves gates pay for the second finish (Clifford+T and the
adder are a few milliseconds slower, QCC's median stays below Qiskit's:
8 ms against 12).

**Learned.**
- Moving a blocker out of the whole runs around the two gates finds
  little: a √X earlier in the run stops it. Cutting the runs where the
  blocker stops, and carrying the gates on its third qubit that only sit
  between in list order, matched pytket on the adder.
- On the adder gathering gives 37 / 85 / 181 CNOTs but 108 / 252 / 540
  gates in all (pytket's own rewrite: 37 / 115 at 8 qubits). The elided
  SWAPs leave generic single-qubit gates; the CX trades did not win them
  back, so the no-more-of-either rule declines it there.
- The CX trades before routing cost the 8-qubit Clifford+T circuit on the
  chip with CZ 100 more gates (phase folding lost its Z rotations); keeping
  diagonal runs diagonal was not enough, running them last was.

**Next.**
1. Win back the adder's gathered CNOTs without the single-qubit cost:
   choose the outer gates of SWAP·U syntheses with their neighbours in view
   (the elided SWAP's blocks are generic where the rest is Clifford+T).
2. Quantum Volume speed (still several times slower than Qiskit at 32
   qubits in single runs).
3. The routing variance the four worse chip rows show: route the logical
   circuit with and without gathering, as with the Pauli network.

## 2026-10-29 — Pauli networks with fewer single-qubit gates, and seeded trials

**What.** `core/src/pauli_network.rs` and `core/src/passes.rs`:

- The local Cliffords that set up each CX of the network are no longer
  undone after it. The rotations still to come are read in the new frame,
  which changes none of their weights, and the frame is paid for once, by
  the final Clifford.
- Each controlled-Pauli gate is placed whichever way round needs fewer local
  Cliffords (X on `a` controlling Z on `b` is a bare CX(b, a)).
- Basis changes use √X (native on IBM-style devices) instead of H: Y to Z is
  one √X, X to Z is S then √X.
- Eight syntheses run side by side, seed 0 deterministic and seeds 1–7
  adding up to 0.05·seed of noise to each candidate's score, so that they
  break the greedy choice's near-ties differently; the cheapest after the
  pipeline is kept (as pytket's GreedyPauliSimp does with `seed` and
  `trials`). The result does not depend on the threads.
- Scoring a candidate gate applies it to the 16 pairs of letters on its two
  qubits once, instead of to every Pauli it touches: about twice as fast,
  which pays for most of the extra trials.
- The network is kept only when it needs no more two-qubit gates and no
  more gates in all than the direct compilation (it used to be kept when it
  won with a CX counted as ten gates).

**Why.** The Pauli network's single-qubit overhead was a known weakness: on
the 16-qubit Hamiltonian it saved CNOTs but left 1285 gates in all against
Qiskit's 972, and at 32 qubits it never won. Undoing each CX's basis change
nearly doubled the network's single-qubit gates, and the greedy synthesis
swung by ±10% in CNOTs with tie-breaking alone. These are fixes of our own
synthesis, verifiable end to end today, rather than a new method.

**Measured** (4-vCPU Xeon VM, `benchmarks/compilation/run.py` and
`device.py`, complete runs). Hamiltonian, two-qubit / all gates:
8 qubits 108 / 451 → 107 / 296 (Qiskit 191 / 474), 16 qubits 325 / 1285 →
313 / 792 (Qiskit 393 / 972, pytket 390 / 887); 32 unchanged at 740 / 1568.
On the chip with CX, 8 qubits 352 / 796 → 349 / 715; with CZ 353 / 1419 →
337 / 1339. No other circuit changed. Geometric means over Qiskit: 0.940×
two-qubit gates (was 0.942×), 0.882× all gates (was 0.916×). On 15 other
random Hamiltonians (6–24 qubits, not used for tuning), against 0.4.0:
0.934× the CNOTs and 0.779× the gates, worse on neither in any. Time on the
Hamiltonian, same machine: 34 / 150 / 406 ms → 36 / 169 / 449 ms at
8 / 16 / 32 qubits.

**Learned.** Tracking each wire's pending local Clifford to score how a
basis change merges with it did not beat counting local gates. Larger
noise (≥ 0.8) loses to the deterministic run. At 32 qubits the network
finds 663 CNOTs against 740, for 1598 gates in all against 1568: kept out
by the no-more-of-either rule.

**Next.**
1. Fewer single-qubit gates still in the network (the 32-qubit Hamiltonian
   is 30 gates away from being kept), e.g. by choosing the final Clifford's
   frame per qubit, or synthesizing commuting sets diagonalized at once
   (phase-polynomial synthesis of each set, as in pytket's PauliSimp with
   sets).
2. Quantum Volume speed (1.8–5× slower than Qiskit on this 4-core VM):
   cheaper two-qubit block resynthesis.
3. The few extra two-qubit gates against Qiskit on the chip (Quantum
   Volume 8, QAOA 16, QFT 16).
