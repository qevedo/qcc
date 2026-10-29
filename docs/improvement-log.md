# Improvement log

One entry per change made by the daily improvement routine: what changed,
why, what it measured, and what to try next. Numbers come from the scripts
in `benchmarks/`.

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
