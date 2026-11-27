# Improvement log

One entry per change made by the daily improvement routine: what changed,
why, what it measured, and what to try next. Numbers come from the scripts
in `benchmarks/`.

See [experiments.md](experiments.md) for the short index: what worked, what
didn't (and why), what was considered, and the current next steps. Read it
before choosing a change, and update it in the same commit.

## 2026-11-27 — Placements chosen for the chip network's rotations (0.4.17)

**What.** `core/src/pauli_network.rs` gains `placement`: a layout of a
circuit's logical qubits on the chip, chosen for its Pauli rotations. It
lowers the sum, over the rotations (as `extract` reads them), of each
one's spanning-tree cost once placed (`tree_cost`, what `synthesize_on`
scores), by first-improvement exchanges of two physical qubits (an
unused one too) until none lowers it. It starts from a given layout, or
from a greedy growth (`greedy_placement`): the logical qubit in most
rotations on the most central physical qubit, then, among the free
qubits next to placed ones, the unplaced qubit sharing most rotations
with that qubit's placed neighbours, the placement of Vuorenkoski and
Meijer-van de Griend ("Shortest Path in Pauli Forest — An Algorithm for
Decomposing Pauli Exponentials to Quantum Circuits", arXiv:2605.03545).
`networks_on_chip` (`core/src/device.rs`) now runs its trials from
SABRE's layouts as before, plus each of them improved and the greedy
one; it keeps the best trial of each ending from SABRE's layouts and
from the new ones apart, and all of them are optimized beside the other
candidates, so the candidates of 0.4.16 are all still there and no
result can be costlier.

Tests: placements are injective layouts of the right length on a line, a
2 × 3 grid and a 7-qubit tree, never costlier than their start and
cheaper on most (not vacuous); a shuffled chain of ZZ rotations on a
line is placed along the chain, every rotation's tree costing 2, from
the greedy start and from the identity. Mutation caught: refusing every
exchange fails both.

**Why.** Next item 1 since 2026-11-23: the chip network's trials only
started from layouts SABRE found for the circuit as written, CX ladders
and all, and the layout was known to matter (first layout only 6562
against 6004 held-out CNOTs). The Pauli Forest paper (read this run)
chooses placement together with architecture-aware Pauli-exponential
synthesis; a placement scored with the synthesis's own tree costs could
be checked end to end today on the held-out set.

**Measured** (4-vCPU Xeon VM; `device.py` full run, all outputs
verified; `check.py` passes; held-out script in the run's scratch space).
- Chip with CX: Hamiltonian 187 / 443 → 186 / 430 (8 qubits), 1063 /
  1894 → 1035 / 1957 (16); with CZ 181 / 730 unchanged, 1030 / 3952 →
  1010 / 3833. Every other row of `device.csv` identical. Geometric means
  over Qiskit level 3: CX 0.879× → 0.877× two-qubit, 0.788× all gates;
  CZ 0.858× → 0.857×, 0.805× → 0.803×.
- Held out (random Hamiltonians built as `circuits.hamiltonian`, seeds
  100n + 7 and 100n + 8, n = 6–16, on the chip): CX 10463 / 21300 →
  10068 / 20220 (−3.8% CNOTs), CZ 10214 / 39643 → 9963 / 38571 (−2.5%);
  none more two-qubit gates; 3 CX rows more gates in all, at most 46
  (n = 15 seed 1507, 876 / 1580 → 825 / 1626). All 44 held-out outputs
  verified by state.
- All-to-all and inside Qiskit unchanged (the plugin does not route).
- Time, paired workers (0.4.16 wheel against this wheel, medians of
  five, alternating order): Hamiltonian-8 189 → 279 ms with CX, 504 →
  577 ms with CZ; -16 1.31 → 2.30 s and 2.89 → 4.31 s (Qiskit about
  0.1 s); other rows 0.88–1.08× (noise on unchanged paths); geometric
  means over the 16 chip rows 1.04× CX, 1.03× CZ.
- Core 125 tests (fmt and clippy clean on Rust 1.100), Python 219
  passed; wasm builds.

**Variants dropped** (details in experiments.md). New placements in place
of SABRE's layouts (greedy alone 10936; greedy and improved 10265, worse
at 8 and 16 qubits); fewer trials on the new layouts (every 2nd 10142,
every 3rd 10197); one best of all layouts per ending (a held-out CZ row
258 → 259 CZs).

**Learned.** SABRE's layouts and the rotation placements find different
good networks; neither alone is best, as with the trial portfolios of
2026-11-07. A pool ranked by an estimate can lose a row when it grows:
rank the new group apart.

**Next.** See experiments.md: chip network time (now the main cost);
a better placement objective or search; the tail's row order.

## 2026-11-25 — The chip network's final Clifford, synthesized on the chip (0.4.16)

**What.** `core/src/pauli_network.rs` gains `State::clifford_on`: on a
device, the Clifford a network leaves at the end is synthesized with
every CX on a coupled pair. Row by row of its tableau, the row's two
images are gathered onto one free qubit `t` (`pair_on`, `pair_from`):
the first along a Steiner tree of the free qubits rooted at `t`
(`steiner_tree`: the shortest path to the closest terminal, again and
again, the Takahashi–Matsuyama heuristic; `gather` clears a letter whose
parent has one with one CX, moves one into an empty parent with two),
made X (or Z); the second, which then anticommutes with it on `t`, along
a tree of the free qubits without `t` onto the neighbour of `t` needing
fewest CXs, cleared by one CX that keeps `t`'s letter. Both orders (X
image first, Z image first) are priced and the cheaper kept. Of every row
and every free `t` that is not a cut vertex of the free qubits (so the
rows left stay on one connected piece of the chip), the one needing
fewest CXs, then fewest gates, goes first; signs are fixed at the end as
before (`fix_signs`, shared with `clifford`). This follows the
architecture-aware synthesis of Clifford tableaus of Winderl et al.
("Architecture-Aware Synthesis of Stabilizer Circuits from Clifford
Tableaus", arXiv:2309.08972). `synthesize_on` now returns both endings
of a trial (the Clifford on the chip, and as before synthesized off the
device), sharing all the rotations before; `networks_on_chip`
(`core/src/device.rs`) cleans and routes both, and passes the best trial
of each ending to the physical pipeline beside the other candidates, so
the candidates of 0.4.15 are all still there and no result can be
costlier.

Tests: the chip-network test checks both endings on a line, a 2 × 3 grid
and a 7-qubit tree whose inner qubits are all cut vertices, with Clifford
circuits on any pairs among the Pauli exponentials, every CX of the
on-chip ending on a coupled pair and every output equivalent; random
Clifford circuits alone on the same chips, all CXs coupled, equivalent,
and fewer CXs than the off-device Clifford routed on more than 30 of 60;
Steiner trees through allowed qubits only, parents before children, and
None when cut off; `connected`. Mutations caught: dropping the cut-vertex
rule; dropping the second CX of a move into an empty parent.

**Why.** Next item 1, measured first: with the final Clifford dropped
altogether (a throwaway switch, not equivalent), 22 random Hamiltonians
of 6–16 qubits on the chip needed 9601 CNOTs instead of 10874 (−11.7%),
Hamiltonian-16 1008 instead of 1078: the routed Clifford tail was about a
tenth of the chip network's cost, more than any other single item on the
Next list could show with a measurement today. Architecture-aware
Clifford synthesis is published and exact, and could be verified end to
end in a run.

**Measured** (4-vCPU Xeon VM; `device.py` full run, all outputs
verified; `check.py` passes; held-out script in the run's scratch space).
- Chip with CX: Hamiltonian 189 / 465 → 187 / 443 (8 qubits), 1078 /
  1903 → 1063 / 1894 (16); with CZ 184 / 742 → 181 / 730 and 1095 / 4198
  → 1030 / 3952. Every other row of `device.csv` identical. Geometric
  means over Qiskit level 3: CX 0.880× → 0.879× two-qubit, 0.790× →
  0.788× all gates; CZ 0.862× → 0.858×, 0.809× → 0.805×.
- Held out (random Hamiltonians built as `circuits.hamiltonian`, seeds
  100n + 7 and 100n + 8, n = 6–16, on the chip): CX 10874 / 21675 →
  10463 / 21300 (−3.8% CNOTs), CZ 10589 / 40489 → 10214 / 39643 (−3.5%);
  none more two-qubit gates (21 of 22 fewer with CX, 19 with CZ); 5 and 5
  with more gates in all, at most 56 (n = 14 seed 1408 with CX, 706 / 1373
  → 677 / 1429). All 48 outputs (with the two benchmark Hamiltonians)
  verified by state.
- All-to-all and inside Qiskit unchanged (the plugin does not route;
  `check.py` passes).
- Time, paired workers (0.4.15 wheel from PyPI against this wheel,
  medians of five): Hamiltonian-8 170 → 235 ms with CX, 476 → 640 ms with
  CZ; -16 1.45 → 1.67 s and 3.60 → 3.96 s (Qiskit about 0.1 s); other
  rows 0.84–1.14× (noise on unchanged paths), geometric means 1.07× CX,
  1.03× CZ.
- Core 123 tests (fmt and clippy clean on Rust 1.100), Python 219
  passed; wasm builds.

**Variants dropped** (details in experiments.md). The chip Clifford
alone, X image first (held-out 10488, but Hamiltonian-8 189 → 191); both
endings ranked together with one optimized (10488 again: the chip ending
wins the ranking nearly always, not the optimized result); the Steiner
search weighted by what `gather` pays (1 per terminal, 2 per other qubit:
identical counts, 1.2× slower).

**Learned.** Dropping a stage to measure its bound took one run and
picked the item; the gap left to it (10463 against 9601) says the tail's
greedy order is now the limit, not the trees. `debug_assert!` does not
run under `cargo test --release`: run new tableau code's tests once in
debug.

**Next.** See experiments.md: a layout chosen for the network itself; the
tail's row order with a lookahead; chip network time.

## 2026-11-23 — Pauli networks synthesized on the chip's pairs (0.4.15)

**What.** `core/src/pauli_network.rs` gains `synthesize_on`: the same
greedy network synthesis on the physical qubits of a coupling graph. Every
CX that sets up a rotation is on a coupled pair; candidates are every
coupled pair a first-layer rotation touches (also with a qubit outside its
support, so a CX can grow a rotation towards its other qubits); a Pauli
costs `tree_cost`, one more than the minimum spanning tree of its support
under the chip's distances (its weight when connected), instead of its
weight, with each Pauli's four possible costs per pair kept in a map by
support; the Clifford left counts 0.2 (`DEVICE_TABLEAU_WEIGHT`, 0.1
off a chip); after 32 gates without a rotation (`DEVICE_PATIENCE`, 4 per
qubit off a chip) the first rotation is reduced by force, its two
closest qubits by a CX if coupled, else one letter moved a step along a
shortest path (two CXs). The final Clifford is synthesized as before and
routed. `core/src/device.rs`: `networks_on_chip` runs the 12 usual network
trials from every distinct initial layout of the routed candidates (48 on
the benchmark Hamiltonians), cleans each (rotation merging, phase folding,
commutative cancellation), routes it from the identity layout with SABRE,
and passes the one with fewest CXs plus SWAP cost to the physical
pipeline beside the other routed candidates; the cheapest is kept as
before. `compile_with(.., on_chip)` lets the tests compare with and
without. `pauli_exponentials` moved to `testing.rs` for the three test
modules.

Tests: `tree_cost` on a line and a grid; 24 networks on a line and a
2 × 3 grid keep their unitary (up to the relabelling) with every CX before
the last rotation on a coupled pair (a mutation allowing all pairs is
caught); on an 8-qubit heavy-hex cell with tails, 8 Pauli-exponential
circuits (CX and CZ) are never costlier with the new candidate, at least 4
cheaper, all on coupled pairs and equivalent up to the layouts (a mutation
in the final-layout composition is caught).

**Why.** The largest measured weakness: the chip's Hamiltonians needed 3–4×
their all-to-all CNOTs (16 qubits: 313 all-to-all, 1321 on the chip),
because the network was synthesized as if all qubits were coupled and
routed after, so the direct compilation won there. Architecture-aware
synthesis of Pauli rotations is published (Paulihedral, Li et al., ASPLOS
2022, arXiv:2109.03371; Tetris, Jin et al., ISCA 2024, arXiv:2309.01905,
both along trees on the device); doing it inside QCC's existing network
greedy keeps CX sharing between rotations and was verifiable end to end
today, while the Next list's items 1–2 (incremental CZ rounds) are larger.

**Measured** (4-vCPU Xeon VM; `device.py` full run, `check.py`).
- Chip with CX: Hamiltonian 335 / 634 → 189 / 465 (8 qubits), 1321 / 2089
  → 1078 / 1903 (16); with CZ 307 / 1099 → 184 / 742 and 1231 / 4201 →
  1095 / 4198. Every other row of `device.csv` identical. Geometric means
  over Qiskit level 3: CX 0.924× → 0.880× two-qubit, 0.811× → 0.790× all
  gates; CZ 0.897× → 0.862×, 0.829× → 0.809×. All outputs verified.
- Held out (random Hamiltonians as `circuits.hamiltonian`, seeds 100n+7
  and 100n+8, n = 6–16, chip, CX and CZ): all 24 fewer two-qubit gates
  (15–35%); CX 8145 / 13551 → 6004 / 12032, CZ 7471 / 26601 → 5799 /
  22044; one row 14 gates more in all (CX n = 16 seed 7: 1273 / 1943 →
  1015 / 1957). All verified.
- All-to-all and inside Qiskit unchanged (`check.py` passes; the plugin
  does not route).
- Time, paired workers (0.4.14 wheel against this wheel, medians of 5):
  Hamiltonian-16 0.49 → 1.39 s with CX, 2.05 → 3.00 s with CZ; -8 1.20×
  and 1.16×; geometric mean over the chip rows 1.11× CX, 1.02× CZ, the
  rest noise (0.68–1.28× on rows whose code path did not change).
- Core 121 tests (fmt, clippy clean on Rust 1.100), Python 219; wasm builds.

**Variants dropped** (details in experiments.md). Weight scoring with
pairs restricted to coupled ones (3700+ CXs routed against 1412 for the
direct one); spanning-tree scoring without neighbour pairs (1295–1321 at
16 qubits); Clifford weight 0.1 / 0.15 / 0.25 / 0.3 / 0.5 / 1.0; patience
4–16, 38, 48, 76; ranking trials on raw CX counts (the cleanup removes
two thirds); trials from the first one or two layouts only.

**Learned.** The greedy on a chip leaves about 3000 raw CXs for 1100 kept:
rank its trials after a cheap cleanup, or the one optimized is a poor
pick. Statevector checks of held-out Hamiltonians must use the decomposed
(Trotterized) circuit, not the `PauliEvolutionGate`, whose `Statevector` is
the exact exponential (all failed, and slowly, until fixed). Dense
unitaries of 12 qubits in a unit test take minutes; 8 is the limit.

**Next.** See experiments.md: a Clifford tail synthesized on the chip; a
layout chosen for the network itself; the network on a chip for the
logical-pipeline candidates of the plugin; CZ time.

## 2026-11-21 — Rotation pairs around two-qubit gates (0.4.14)

**What.** `core/src/pauli_pairs.rs` gains `insert_rotation_pairs`, run once
on every final result (`with_rotation_pairs` in `passes.rs`, at the end of
`until_fixed_point`, which all-to-all, plugin and each chip candidate go
through):
- Where a whole axis commutes with a two-qubit gate on one qubit (tested
  on the gate's matrix with a rotation by a generic angle: Z on a CX's
  control and on a CZ, X on a CX's target), any rotation D about it gives
  G = D† G D. The wire DP of Pauli pairs now has, at such a gate, the
  identity, the Pauli and rotations at angles taken from the runs in
  reach: the two nearest non-empty runs on each side that the axis
  reaches through gates of the same axis (`WINDOW`), each offering the
  angle that removes its outer rotation in that axis (from its ZYZ angles
  in that basis), shifted also by π and −π/2 (the other outer angles of
  the {rz, sx} form), or its whole angle if it is a rotation about the
  axis alone. Gates with only some Paulis commuting (ECR) keep those.
- Pauli pairs inside the pipeline are unchanged (same code with rotations
  off; the 234 compilations below are byte-identical to 0.4.13 before the
  last pass). Rotations run only on the final result, so every keep rule
  and candidate choice sees exactly what it saw in 0.4.13; the result is
  kept only if no costlier (a guard; it never triggered).
- Speed: between two Z-axis gates, a run priced with Rz(α) on one side
  and Rz(−β) on the other has ZYZ angles (θ, φ + α, λ − β), and where its
  count is `base + [λ + dλ ≠ 0] + [φ + dφ ≠ 0]` (`Basis1q::split_count`,
  new in `one_qubit.rs`) the best choice before is found once for all
  choices after (K instead of K² prices, same ties). Elsewhere such runs
  are priced from their angles (`count_from_angles`) and the rest by
  matrices. Euler angles once per run, flat choice and DP buffers, and
  predecessors that already cost as much as the best total are skipped.
  QV-32 with CZ: 2.86G → 1.00G instructions over 30 calls; QFT-32 with CX
  (runs bordering CX targets, matrix path) 9.4M per call against 2.7M for
  the Pauli pass. Wall time of the pass on the 32-qubit outputs 0.2–0.8
  ms, Quantum Volume 3.2 ms (CZ) and 5.9 ms (CX).
- Tests: a Z rotation through two CX controls and a CZ into the next run
  (13 → 12 gates; a lone Rz between runs, two fewer), an X rotation hidden
  on both sides of two CX targets is found, 600 random circuits (CX, CZ,
  ECR, {rz, rx, ry}) never costlier than Pauli pairs and cheaper on more
  than 50, the separable pricing gives identical outputs to the generic
  one on 400 circuits (a sign mutation there is caught only by this test),
  and `split_count` equals `count_from_angles` on 1000 samples. The
  pipeline test now asserts that a CX compile is the first pipeline plus
  `with_rotation_pairs`, exactly.

**Why.** Next item 8 (Pauli pairs beyond Paulis). The research turned up
nothing measured to beat QCC's weak spots today (see experiments.md,
Considered); the CZ-time and routing items need much larger work. This one
could be verified end to end in a run and, run last, cannot make any
circuit costlier. Qiskit's `Optimize1qGatesSimpleCommutation` moves gates
one commutation at a time into existing runs; this solves each wire
exactly over a candidate set that contains the Pauli solution.

**Measured** (4-vCPU Xeon VM; `run.py`, `device.py`, `plugin.py`, full runs).
- All-to-all: Clifford+T 274 → 271, 484 → 476, 1009 → 988 (pytket 280,
  495, 1039); QAOA 312 → 310, 615 → 612; Hamiltonian 286 → 284, 786 →
  781, 1504 → 1499. Gates 0.866× → 0.863× of Qiskit, two-qubit gates
  unchanged (0.923×). All outputs verified.
- Chip with CX: adder 151 → 148, 402 → 395, QAOA-16 403 → 401,
  Hamiltonian 636 → 634, 2091 → 2089, Clifford+T 507 → 506, 1207 → 1206;
  with CZ: Clifford+T-16 2668 → 2664, Hamiltonian 1103 → 1099, 4203 →
  4201, adder-16 614 → 612, QV-8 690 → 689. Inside Qiskit: adder 160 →
  157, 405 → 398 (CX), Hamiltonian 1286 → 1283, 4421 → 4415 and QV-16
  4266 → 4259 (CZ). Chip gates against Qiskit 0.813× → 0.811× with CX,
  0.829× with CZ (unchanged at three decimals); inside Qiskit 0.842× →
  0.839× with CX, 0.830× with CZ. All outputs verified.
- Over 234 compilations (the benchmark's 112 and 122 held out: sizes 6–24,
  all-to-all and chip, CX and CZ): 67 with fewer gates (32 / 118 gates on
  the benchmark, 35 / 150 held out), none with more gates, no two-qubit
  count changed.
- Time, paired workers (0.4.13 wheel against this wheel, alternating):
  all-to-all 1.03× in geometric mean (medians of 11); chip 1.00× (medians
  of five). A harness check, old against old, read 1.002×.
- `check.py` passes; core 118 tests (clippy clean on Rust 1.100), Python
  219; wasm builds.

**Variants dropped** (details in experiments.md). Rotations inside the
pipeline (12 rows with more two-qubit gates: the adder lost gathering,
181 → 196 CNOTs at 32 qubits); window 1 (4 rows a gate worse) or 4 (one
gate better, more time); the pass twice (4 gates over 234 compilations);
Z only or X only (each loses most of the other's rows); fewer angle
offsets (14–48 rows worse); a sorted cut-off in the DP (slower).

**Learned.** Run a pass with a strict never-worse guarantee after every
decision, not inside, when its gains would move what keep rules compare.
A fast path that only changes prices passes equivalence and never-worse
tests even when wrong; test it against the slow path directly.

**Next.** CZ device time (incremental rounds) and routing quality stay
first (experiments.md Next); for these pairs, letting rotation gains inform the keep rules
without the fall-back (chip CX Hamiltonian-8 got 333 two-qubit gates that
way, against 335).

## 2026-11-19 — Pauli pairs around two-qubit gates (0.4.13)

**What.** New pass `core/src/pauli_pairs.rs` (`insert_pauli_pairs`), run
in `passes.rs`:
- If a Pauli P on one qubit commutes with a two-qubit gate G (checked on
  G's matrix, up to phase: X on a CX's target, Z on its control, Z on
  either qubit of a CZ), then G = P G P. Putting a P on each side merges
  one into each neighbouring run. Along a wire this carries a Pauli from
  run to run across any number of such gates. Pairs on different wires
  are independent, so each wire is solved exactly: a dynamic programme
  over the wire's chain of runs, four states (I, X, Y, Z) at each gate,
  each run priced by `count_1q` of P_after · R · P_before (as it is when
  both are I). A wire changes only if strictly cheaper; only changed runs
  are resynthesized, and a run is created where a wire had none if that
  pays.
- It runs at the end of `tail` in every pipeline whose result is final
  (`Last::Yes`, `Last::Carrying`), so the gathering keep rule and the
  carrying rule compare circuits that have had it. For CX-like devices
  other than CX it runs instead once after `until_fixed_point`'s rounds.
- Masks of commuting Paulis are cached per gate kind; empty slots use a
  precomputed table (speed only, outputs identical).
- Tests: an X carried through two CX targets to a generic run; no Pauli
  crosses a CX control or a measurement; two Xs cancel and a Z crosses a
  CZ; 450 random native circuits (CX, CZ, ECR) keep their unitary, never
  cost more, keep their two-qubit gates, and more than 50 get cheaper.
  Two mutations (the product's order, a wrong commutation mask) are caught.

**Why.** Next item 7 was the adder's gathered CNOTs (37 / 85 / 181
against 40 / 92 / 196 kept), refused because the gathered circuits had 3,
7 and 15 more gates in all. Counting gate names in the refused circuits
showed the excess was exactly their `x` gates: each an X left by an
elided SWAP, alone after a CX it is the target of. Within a few gates
along the wire there is a run that would absorb it for free, but no pass
looked further than one gate (`shift_through_blocks` trades one lone CX
at a time; Qiskit's `Optimize1qGatesSimpleCommutation` likewise moves one
commutation at a time into existing runs). A Python prototype on QCC's
outputs (20 minutes) showed the same Xs everywhere: the 32-qubit Heisenberg
chain had 66. Exact and cheap, it beat the Next list's other items, which
need a decision on a keep rule (carrying) or large work (incremental rounds).

**Measured** (4-vCPU Xeon VM, `run.py`, `device.py`, `plugin.py`, full runs).
- All-to-all with CX: the adder 40 / 105 → 37 / 104, 92 / 245 → 85 / 244,
  196 / 525 → 181 / 524 (gathering now kept); Clifford+T 278 → 274,
  499 → 484 (pytket 495), 1030 → 1009; Heisenberg 273 → 267, 573 → 559,
  1173 → 1143; QAOA 159 → 156, 317 → 312, 634 → 615; Hamiltonian
  790 → 786, 1511 → 1504. Geometric mean against Qiskit 0.932× / 0.874× →
  0.923× / 0.866×.
- Chip with CX 0.816× → 0.813× of Qiskit's gates (Heisenberg-16 571 →
  559, Hamiltonian-16 2099 → 2091); with CZ 0.830× → 0.829× (Clifford+T-16
  2677 → 2668). Inside Qiskit 0.844× → 0.842× with CX, 0.835× → 0.830× with
  CZ (adder-16 613 → 597). No two-qubit count changed off the adder, and
  no row of any table got costlier.
- Held out (other sizes and seeds, 6–24 qubits, all-to-all and chip, CX
  and CZ): 132 compilations, 6 with fewer two-qubit gates (all adders,
  7–8%), 65 with fewer gates, none worse in either.
- Time, paired workers (0.4.12 wheel against this build, alternating,
  medians): 1.04× all-to-all, 1.03× chip CX, 0.98× chip CZ in geometric
  mean; Clifford+T-32 18.5 → 19.4 ms, Hamiltonian-32 131 → 132 ms.
- `check.py` passes; core 113 tests, Python 219; wasm builds.

**Variants dropped.** Inside the CZ rounds (some rows better, e.g. chip
Clifford+T-16 2649 gates, but others reached worse fixed points: chip
Hamiltonian-16 1226 / 4229 against 1231 / 4208, all-to-all QAOA-8 223
against 213). In every pipeline, also before routing: chip CX adder-16
146 / 402 → 151 / 456, Hamiltonian-16 1321 → 1339 CNOTs.

**Learned.** A refused keep-rule trade can hide a missing pass: count
the gate names of the refused circuit. A pass that never makes its input
costlier can still make the CZ rounds end worse; such passes go after
`until_fixed_point`. A prototype that prices unchanged runs differently
from the real pass (here through a 1-qubit compile) overestimates.

**Next.**
1. Incremental CZ rounds (unchanged; CZ Hamiltonian-16 on the chip 2.0 s).
2. Pauli pairs inside the CZ rounds, with both trajectories compared.
3. Carrying beyond the strict rule: remeasure with Pauli pairs.

## 2026-11-17 — One map update per gate in rotation merging (0.4.12)

**What.** `core/src/pauli_network.rs`:
- `apply_all_to_map`: the Clifford map that rotation merging and the Pauli
  network's extraction keep (`C† X_q C`, `C† Z_q C`) is updated once per
  gate instead of once per Clifford the gate is read as. Appending
  g₁, …, gₖ maps each generator P to `C† g₁† ⋯ gₖ† P gₖ ⋯ g₁ C`, so P is
  conjugated by gₖ† first and g₁† last, and its image read once. √X was
  three updates (H S H), X four, Y six, CZ and SWAP three, Z rotations by
  π/2, π and 3π/2 one to three; each is now one.
- `read` returns its steps inline (`Step::Cliffords([Gate; 6], n)` or
  `Step::Rotation`) instead of a `Vec` per instruction.
- Test: `one_update_per_gate_is_the_same_map` (every multi-step gate
  `read` gives and S², S³, on two qubits of three in both orders, on maps
  scrambled by random Cliffords, against one update per Clifford). A
  mutation that drops the inverses fails it; the step order alone cannot
  be caught, as every gate's steps equal their reverse up to phase.

**Why.** Next item 1 was the CZ chip time (16-qubit Hamiltonian 2.44 s
against Qiskit's 0.11–0.15 s). The suggested remedy, running the rounds
of `until_fixed_point` only on the best routed candidates, was measured
first and does not hold (see **Didn't work** in `experiments.md`): the
first pass does not predict which candidate ends cheapest. A callgrind
profile of one candidate's rounds (the routed circuit dumped from the
device path and run through `passes::compile` in a throwaway example)
was spread out: two-qubit synthesis 28%, gathering 11%, M4 products 10%,
rotation merging 8.9%, of which 7.0% was `apply_to_map`, mostly on the
√X and CZ gates the CZ pipelines are full of. That part could be removed
exactly.

**Measured** (4-vCPU Xeon VM).
- Outputs: 220 compilations byte-identical to 0.4.11 (every circuit
  all-to-all and through the Qiskit plugin with CX and CZ; the chip with
  CX and CZ including 27 held-out circuits per gate set; 31 held-out
  circuits on the chip with CZ and 23 on a 4×5 grid with CX, other seeds).
- Callgrind, one CZ candidate of the 16-qubit Hamiltonian on the chip:
  15.30 → 14.74 G instructions (−3.7%); rotation merging 1.37 → 0.81 G.
- Paired timings (two long-lived workers, 0.4.11 wheel and this build,
  alternating, medians of 9 all-to-all and 5 on the chip): all-to-all
  1.00× with CX and 0.974× with CZ in geometric mean; chip 1.001× with CX
  and 0.983× with CZ, the CZ Hamiltonian-16 2.77 → 2.57 s (0.93×). Single
  circuits of a few milliseconds vary by ±15% between runs here.
- `check.py` passes; core 109 tests, Python 219.

**Learned.**
- Pruning routed candidates on their first-pass cost does not work (97
  CZ compilations, 61 with several candidates): no margin up to 7% keeps
  every result, and at 10%, the first that does, only 10% of the rounds
  go (145 of 1518). Pruning in lockstep after every round is no better
  (15% at 10%). The winner of
  the 14-qubit held-out Hamiltonian started 2.1% behind another
  candidate.
- In the rounds neither pipeline can be dropped: the direct one won 368
  rounds (365 strictly), the one through CX 139.
- `cx_count` in gathering is a full KAK per call (5.3% of a CZ
  candidate's instructions); a memo by the unitary's bits would hit 44%
  of calls, at most 2.3%, not worth a cache passed through gathering.

**Next.**
1. The CZ rounds' time is spread over the whole pipeline; what remains
   large is resynthesis (synthesis misses 28%, `dress` and `kak`) and
   gathering's KAKs. Incremental rounds (re-run only the parts of the
   circuit the last round changed) would cut it by much more than any
   one pass, but need a way to bound where a change can reach.
2. Carrying beyond the strict rule (unchanged).
3. Quantum Volume speed all-to-all (QV-32: 96 ms, the core 78 of it).

## 2026-11-15 — Carrying gates between blocks on the same pair (0.4.11)

**What.**
- `core/src/commuting_blocks.rs`: `carry_through`. For each two-qubit gate,
  the next two-qubit gate on its pair is found, with at most three
  blockers (two-qubit gates on one of the pair's qubits and a third)
  and any single-qubit gates on the pair between. The end of the first
  block, from one of its last three two-qubit gates, is carried forward:
  it crosses each single-qubit gate g as g D g† and each blocker only if
  it commutes with D as it stands there (an exact 8×8 check), and is put
  down before the next block as its own gates between the inverse of the
  single-qubit gates it crossed and those gates again (R† D R), for the
  resynthesis to merge. Likewise the start of the second block, carried
  back. A move is made when the two blocks then need fewer CX-like gates
  (one KAK each, only for parts that cross).
- `core/src/passes.rs`: `pipeline` takes `last: Last` (`No`, `Yes`,
  `Carrying`) instead of a bool. With `Carrying`, `finish` carries after
  block resynthesis, resynthesizes again and finishes both circuits,
  keeping the carried one only if it needs no more two-qubit gates and
  no more gates in all. The choice between gathering and not is now made
  on the costs without carrying (`plain_cost`). For CX-like non-CX
  devices, `until_fixed_point` runs its rounds without carrying until
  they stop shrinking, then with carrying (skipped when nothing in the
  result can be carried); for CX the first final pipeline carries.
  Routing's logical pipeline (`Last::No`) never carries.
- Tests: `a_cx_carried_through_a_hadamard_joins_the_next_block`, its mirror
  backwards, `nothing_is_carried_through_a_gate_it_does_not_commute_with`,
  `carrying_keeps_random_circuits_equivalent` (300 random Clifford+T-like
  circuits, unitaries checked) and, end to end,
  `a_cx_is_carried_to_the_next_block_on_its_pair` (3 CXs → 2 with CX; for
  CZ the rounds already found it).

**Why.** Next item 4 of the last run: a block on {a, b}, a gate on {b, c},
a block on {a, b} again, where part of the first block commutes with the
middle gate once moved through the single-qubit gates between
(CX(15,1) H(1) CX(1,7) CX(15,1) is two CXs, not three). Gathering moves
only whole gates, unchanged. A Python prototype applied the move to
QCC's final outputs and recompiled them: all-to-all it saved 2–6 CNOTs on
QAOA, Clifford+T and the Hamiltonian (with more single-qubit gates on
some); on the chip, where every gate stays on its coupled pair, 29 CNOTs
on Clifford+T-16, 10 on the 8-qubit Hamiltonian, 14 on the 16-qubit one
(CX), 23 with CZ. This beat the other Next items as a measured quality gap
that is coupling-safe, exact and testable today; the speed items (CZ
rounds, Quantum Volume) make no circuit smaller. Searched beforehand:
`ThreeQubitSquash` and phase-polynomial resynthesis were dropped on
2026-11-13; ZX-calculus extraction is not coupling-aware.

**Measured** (4-vCPU Xeon VM; `run.py`, `device.py` and `plugin.py`
rerun in full, every output verified).
- All-to-all: Clifford+T 16 150 / 505 → 148 / 499 (pytket 148 / 495),
  32 310 / 1034 → 308 / 1030 (pytket 310 / 1039); geometric mean over
  Qiskit 0.933× / 0.875× → 0.932× / 0.874×. Every other row unchanged;
  `check.py` passes.
- Chip with CX: QAOA-8 94 / 207 → 92 / 207, adder-16 148 / 405 →
  146 / 402, Hamiltonian-8 345 / 648 → 335 / 639, -16 1336 / 2100 →
  1321 / 2099 (geometric mean 0.928× / 0.818× → 0.924× / 0.816×). With
  CZ: Clifford+T-16 692 / 2677 → 691 / 2677, Hamiltonian-8 310 / 1119 →
  307 / 1103, -16 1256 / 4454 → 1231 / 4208 (0.899× / 0.834× → 0.897× /
  0.830×). Inside Qiskit: adder-16 CX 159 / 409 → 157 / 405; CZ
  Clifford+T-8 242 / 848 → 240 / 840, -16 741 / 2455 → 735 / 2439,
  Hamiltonian-8 388 / 1318 → 354 / 1288, -16 1312 / 4416 → 1294 / 4434
  (18 fewer CZs, 18 more gates in all: the CZ rounds keep a round with
  fewer two-qubit gates). No row has more two-qubit gates.
- Held out: 82 compilations (Hamiltonians of 6–24 qubits, Clifford+T of
  6–24, QAOA 6–18, adders 10–14, QFT 10–12, Heisenberg 10–12;
  all-to-all and on the chip up to 18 qubits, CX and CZ) against 0.4.10:
  fewer two-qubit gates on 22, more on none (124 fewer, 415 fewer gates
  in all); two CZ chip rows have 4–6 fewer CZs and 2–12 more gates in
  all.
- Time, against 0.4.10 side by side, interleaved (medians of 9
  all-to-all, 5 on the chip): 1.05× all-to-all and 1.07× on the chip in
  geometric mean. Where moves are found it costs a second resynthesis
  and finish: Clifford+T-32 all-to-all 13.7 → 26.8 ms, the adder 1.14–
  1.21×, the Hamiltonian 1.15–1.19×; the CZ chip Hamiltonian-16 1.46 →
  2.44 s, where the carrying rounds keep removing CZs (1256 → 1241 →
  1235 → 1231 → 1231 / 4208 over four rounds on the winning routing).
  Before the two speed-ups (KAKs only for parts that cross; no carrying
  rounds when nothing can be carried) it was 1.19× all-to-all, with every
  CZ circuit paying one more round.

**Learned.**
- Placement mattered more than the move. Carrying before routing, or in
  every pipeline, changed what SABRE and the CZ rounds saw and made some
  rows worse; on final results only, after resynthesis, it never did.
- The gathering keep rule ("gathered if it dominates, else not") was not
  monotone under a pass that improves both sides: an all-to-all
  Clifford+T-16 result got 2 more CXs until gathering was judged on the
  costs without carrying.
- The CZ rounds are trajectory-dependent: carrying from the first round
  ended Clifford+T-16 at 695 CZs instead of 692; plain rounds first, then
  carrying rounds, is never worse.
- The biggest wins trade single-qubit gates for CXs (chip CX
  Clifford+T-16: 704 / 1251 against 733 / 1209) and are refused by the
  no-more-of-either rule; judging each move by synthesized gate counts
  did not recover them.

**Next.** Chip CZ time (rounds on every routed candidate), carrying
within the rule with fewer single-qubit gates, Quantum Volume speed.

## 2026-11-13 — Pauli-network scoring from tables, shared KAKs in symmetric blocks (0.4.10)

**What.**
- `core/src/pauli_network.rs`: `candidate_tables()` works out once what
  each of the nine candidate controlled-Pauli gates (σ controlling τ) does
  to the 16 pairs of letters on its two qubits, and how many local
  Cliffords it needs. `best_gate` used to build each candidate's gates
  (three small `Vec`s and a concatenation) and apply them to the 16 pairs
  for every pair of qubits and every gate placed. `State::scored` reads
  each Pauli's weight and discount once per step instead of once per pair,
  and `State::changes` sums the nine candidates' scores in one pass over
  those Paulis, each sum in the old order and from −0.0 as
  `Iterator::sum` does, so every score is the same float. Only the
  winner's gates are built. `Gate::Sx` conjugates in one step (the X bit
  takes the Z bit, the phase gains three Z bits) instead of as H, S, H.
- `core/src/symmetric_blocks.rs`: each pass keeps the KAK decomposition of
  every block core it has seen, by the matrix's bits (`key_2q`, now
  `pub(crate)` in `lower.rs`); and for a lone CX-like block, a side whose
  two runs already cost nothing skips the search for a better angle, which
  could not beat it.
- New tests: `sx_is_h_s_h` (every letter and phase against the three
  gates), `candidate_tables_match_the_gates` (the tables against the gates
  placed on qubits 5 and 2 next to a third letter).

**Why.** Next items 2 and 3 were compile time on Hamiltonians: 0.33 s at
32 qubits against Qiskit's 17 ms. Before choosing, two quality leads were
measured and dropped (see experiments.md): phase-polynomial resynthesis
of CNOT+Rz regions (GraySynth + PMH) saved no CNOT on any benchmark
output, and three-qubit resynthesis (pytket's `ThreeQubitSquash` on QCC's
output) saves 2–3 CNOTs all-to-all and its chip "savings" put most of
their CXs on uncoupled pairs. A callgrind profile of the 32-qubit
Hamiltonian in a Rust driver (14.4 billion instructions) put 83% in the
Pauli network's greedy synthesis, 45% in `Gate::conjugate` called from the
candidate scoring and 20% in allocation, all of it exact work that could
be done once.

**Measured** (4-vCPU Xeon VM). Outputs: 142 compilations (every benchmark
circuit all-to-all with CX and CZ, on the chip with QCC alone and inside
Qiskit, plus 15 held-out random Hamiltonians of 6–24 qubits with CX and
CZ) byte-identical to 0.4.9's; `run.py` and `device.py` rerun in full
give the same counts on all 96 + 96 rows. Instructions (callgrind, core
alone, 32-qubit Hamiltonian, billions), step by step: tables and the
one-step √X 14.4 → 4.9, nine sums in one pass 4.9 → 4.1, weights read
once per step 4.1 → 3.8, KAK memo 3.8 → 3.5, skipping free sides
3.5 → 3.4; the network synthesis itself 12.0 → 1.4.
Time against 0.4.9, side by side, best of two medians of 11, all-to-all:
Hamiltonian 8 / 16 / 32 with CX 31.4 → 18.5, 137.2 → 51.0, 340.2 →
110.2 ms (0.59×, 0.37×, 0.32×); with CZ 51.4 → 36.4, 175.7 → 90.8,
444.8 → 188.8 ms (0.71×, 0.52×, 0.42×); the other 42 rows 0.94× in
geometric mean (three that read 1.12–1.16× were 0.88–0.97× in medians of
101). On the chip (best of two medians of five): 0.96× with CX and with
CZ; Hamiltonian 16 438 → 376 ms with CX, 1.45 → 1.37 s with CZ, where
SABRE (24% of the instructions), the repeated CZ rounds (189 pipelines)
and two-qubit synthesis (20%) dominate.

**Learned.** Profile before choosing among speed items: the Next list
guessed "stop a trial early" and "cache templates' KAK", but the cost
was the scoring loop. Scores are compared with a 1e-9 margin, but keeping
each float sum in its old order made the outputs provably identical
rather than merely very likely so. On the CZ chip path no single pass
dominates: 189 pipeline runs share the time between resynthesis (11%),
moving gates through blocks (7%), single-qubit merging (4%) and rotation
merging (4%), and SABRE takes another 24%.

**Next.** Chip CZ time (pipeline rounds and SABRE), the pyo3 crossing on
Quantum Volume, routing quality (LightSABRE's relative scoring), the
cross-pair rewrite found in the 3-qubit analysis (below in experiments.md).

## 2026-11-11 — Instructions without heap memory, a faster hasher, and no GC storms (0.4.9)

**What.**
- `core/src/ir.rs`: an `Instruction`'s name is a `Name` (a
  `Cow<'static, str>` borrowed from a list of the 50 gate names the
  compiler knows, owned only for opaque gates), and its qubits and
  parameters are `Inline<T, 3>`: up to three items stored inline, a `Vec`
  beyond. Both read, compare, hash and print as the `str` and slice they
  hold, so almost every pass compiled unchanged; `Instruction::gate` takes
  anything that converts (arrays, slices, `Vec`s), and its callers now pass
  arrays. `one_qubit::Op` and the two-qubit synthesizer's gate keys use the
  same types.
- `core/src/hashing.rs`: `FastMap`, a `HashMap` with FxHash (rustc's
  multiply-and-rotate hash), for the internal maps of phase folding,
  two-qubit synthesis, the passes, lowering's memos and SABRE's `Savings`.
- `python/qevedo/compiler/passes/native.py`: `paused_gc()` disables
  Python's cyclic collector while instructions cross to and from the core
  (also around the Qiskit plugin's crossing), and restores it after.

**Why.** Quantum Volume at 32 qubits took 243 ms in the last full run
against Qiskit's 49 ms; Next item 1 named allocation as the main serial
cost. A callgrind profile of the core alone (1.09 billion instructions)
confirmed it: malloc, free and friends were 27% of the instructions, and
two thirds of the allocations came from copying, relabelling and building
instructions (`Instruction::clone` 150 000 times, `relabel` 130 000,
`one_qubit::rot` 107 000). SipHash was another 3.5%. Timing the Python
side piece by piece found that building the tuples and `Instruction`s
cost only about 5 ms each way, but the native pass took 95–103 ms against
75 ms for the bare core call: with the collector disabled it took 75–80.
CPython starts a collection every 700 new objects, and the crossing makes
about 45 000, so some of those collections walk the whole heap (Qiskit's
modules included). Qiskit's own Rust core stores standard gates as an
enum and parameters inline for the same reason. These are exact speed-ups
whose outputs can be compared byte for byte, of the biggest measured
weakness against Qiskit; the quality candidates (three-qubit
resynthesis, routing search) promised less certain gains.

**Measured** (4-vCPU Xeon VM).
- Outputs: 160 compilations (24 circuits all-to-all with CX and CZ, 16
  on the heavy-hex chip with CX and CZ, each also through the Qiskit
  plugin) identical byte for byte to 0.4.8's, including output
  permutations; `run.py`, `device.py` and `plugin.py` rerun in full, every
  gate count unchanged, every output verified; `check.py` passes.
- Core alone, QV-32: 1.09 → 0.91 billion instructions (callgrind);
  allocation 27% → 10% of them; 93 → 59 ms (median of 11, Rust binary).
- End to end against 0.4.8 side by side (best of two medians of 11,
  all-to-all, CX): 0.85× in geometric mean over 24 circuits; QV 8/16/32
  6.9 → 5.6, 30.0 → 22.7, 149 → 92 ms; Clifford+T 32 16.6 → 12.9 ms;
  Heisenberg 32 13.3 → 10.4 ms; Hamiltonian 32 358 → 326 ms; GHZ (0.3–0.6 ms)
  0.96–1.04×, noise.
- Chip (best of two medians of five): 0.93× with CX, 0.83× with CZ; CZ
  Hamiltonian-16 1.79 → 1.42 s, CZ QV-16 424 → 335 ms.
- Full benchmark run: QCC's median 6 ms against Qiskit's 10; faster than
  Qiskit on 20 of 24 circuits (was 18); QV-16 20 against 22 ms, QV-32 89
  against 50 ms (was 243 against 49). Chip medians 24 ms with CX (Qiskit
  31) and 77 ms with CZ (Qiskit 32).

**Learned.**
- The Python wrapper's cost was the garbage collector, not the
  dataclasses: earlier runs (2026-11-01, -03, -06) blamed the 9,000
  `Instruction`s, which cost about 5 ms. Time Python glue with `gc`
  disabled as well as enabled before optimizing it.
- Changing a core type behind `Deref` to the old type's slice or `str`
  kept the diff small: 14 compile errors in the library after the type
  change, all mechanical.
- The callgrind numbers include the worker threads' synthesis; the wall
  time gain on QV-32 (0.62×) was larger than the instruction count's
  (0.83×) because the main thread's serial passes were the allocating
  ones.

**Next.**
1. QV-32 is still 1.8× Qiskit. The crossing itself: the core call takes
   70–75 ms from Python against 59 ms in a Rust binary (pyo3 builds a
   list, a tuple and a string per instruction; interned names or a flat
   array would cut it). Then the synthesis: `dress`, `Shape::of`, 4×4
   products, `hypot` and `sincos` are the top of the profile now.
2. CZ device time (Hamiltonian-16 1.42 s against Qiskit's 0.15 s).
3. Hamiltonian compile time (0.33 s at 32 qubits against 17 ms).

## 2026-11-09 — Faster SABRE, and a wider routing search on devices (0.4.8)

**What.** `core/src/routing/sabre.rs` makes the same routing choices at
2.1–2.9× fewer instructions:
- `Dag`: the dependency graph (as flat successor lists), each node's kind
  and each gate's matrix (and its operands-reversed form) are worked out
  once per circuit and shared by every routing of it: a device search
  routes the same circuit hundreds of times.
- The extended set reuses scratch arrays with generation stamps instead of
  copying the predecessor counts and allocating a `seen` array as long as
  the circuit at every SWAP (15% of the time was `memcpy`).
- A routing is a `Trace` of node indices and SWAPs; only the routing kept
  becomes `Instruction`s (before, every routing cloned every name and
  parameter list: malloc/free were a quarter of the time).
- `Savings`: the CNOTs a SWAP saves after a block (two KAK CNOT counts)
  memoized by the block's unitary, bit for bit, across the routings of
  one start; a single-qubit gate is folded into a block with the two
  nonzero terms of `(g ⊗ I)·m` in the full product's order (`left_1q`),
  and a gate the other way round by permuting entries (`reversed`), so
  the numbers are the same as the full products'.
- SWAP scores read the front's and extended set's physical pairs once per
  SWAP, sum distances as integers (exact, so the same floats), and
  `CouplingGraph` keeps distances in one flat array.

`core/src/device.rs`: `routing_options` now has 16 swap trials per layout
(was the default 4, `USUAL_SWAP_TRIALS`); `layout_and_route_and_first`
also returns the routing the first 4 swap trials of each layout would have
given, when it differs, and `compile_for_device` optimizes both and keeps
the cheaper. Tests: `a_single_qubit_gate_after_a_block_is_the_full_product`,
`a_gate_the_other_way_round_is_conjugated_by_swaps`,
`the_savings_memo_changes_nothing`,
`the_first_swap_trials_choose_as_a_narrower_search`,
`the_wide_search_is_never_worse_than_the_usual_one` (it also requires some
circuit to get strictly cheaper), and the device equivalence check with the
device options on undirected CX and directed CZ chips.

**Why.** The ledger listed QV-8 on the chip (111 CNOTs against Qiskit's 108)
as QCC's only row worse than Qiskit. Investigating it showed that row is
seed luck (Qiskit gives 108 or 111 depending on its seed; QCC's SWAP-cost
estimate is exact there and no search width found 108), but also that
SABRE's results vary widely between tie-breakings (Qiskit: 798–840 on
QV-16) and that routing was most of a chip compilation's time with CX
(QV-16: ~160 of ~210 ms). LightSABRE (Zou et al., arXiv:2409.08368) makes
the same point: cheap routing lets you exploit SABRE's randomness with
more trials. A faster router with byte-identical outputs is verifiable
end to end, and spending the time saved on more trials gives gates on
every family that routes.

**Measured** (4-vCPU Xeon VM). With the old options, all 112 benchmark
compilations (all-to-all, chip and inside Qiskit, CX and CZ) are byte
for byte what 0.4.7 gives, and the multicore example's counts are
unchanged. Callgrind on one device search (32 starts, 4 swap trials) of
the 16-qubit circuits: Quantum Volume 4.01 → 1.58 G instructions,
Hamiltonian 3.90 → 1.54, Clifford+T 2.80 → 1.33, QFT 1.94 → 0.66; with 16
swap trials 4.04, 3.47, 3.14 and 1.42 G, about what 4 used to cost.
`device.py` in full: QCC alone 0.934× → 0.928× Qiskit's two-qubit gates
and 0.820× → 0.818× all gates with CX, 0.904× → 0.899× and 0.839× →
0.834× with CZ; seven rows cheaper, none dearer (QV-16 822 → 783 with CX
and CZ, Clifford+T-16 752 / 1251 → 733 / 1209 CX and 708 / 2790 → 692 /
2677 CZ, Hamiltonian-16 1368 / 2106 → 1336 / 2100 CX and 1292 / 4572 →
1256 / 4454 CZ, QAOA-16 CX 180 / 409 → 179 / 403); every output passes
its state check; Qiskit's and Qiskit + QCC's rows are unchanged (Qiskit
routes there). 40 held-out circuits (QV 10–16 qubits seeds 101/202,
Clifford+T, Hamiltonians and QAOA of 10 and 14 qubits seeds 7/8, CX and
CZ): 19 213 → 18 936 two-qubit gates and 66 134 → 64 817 gates; fewer
two-qubit gates on 19, more on none; three have a few more single-qubit
gates for fewer two-qubit ones. Time against 0.4.7 side by side (best of
two medians of five, interleaved): 0.86× with CX and 0.97× with CZ in
geometric mean; slower where both routings go through the CZ rounds:
Hamiltonian-16 CZ 1.70 → 2.33 s, QV-16 CZ 1.22×. All-to-all compilation
does not route, so `run.py`'s results and `check.py` are unchanged.

**Learned.** Wider search alone was not monotone in the final count (16
swap trials alone: 3 held-out circuits worse; 8 or 32: benchmark rows
worse), because the SWAP-cost estimate ignores single-qubit gates, CZ
rounds and network Cliffords; optimizing the old choice beside the new
makes it monotone. A block-aware lookahead, more layout trials and top-k
routings were tried and dropped (experiments.md). An experiment switch
left in `routing_options` made a later edit miss silently; the
never-worse test caught it because it also demands a strict improvement.

**Next.** All-to-all Quantum Volume speed with the same recipe (record,
write out once); CZ device time (the rounds now run twice on some
circuits); routing quality from LightSABRE's relative scoring and block
routing.

## 2026-11-07 — Pauli-network trials judged one by one, and four that charge more for local Cliffords

**What.** `core/src/passes.rs`: the Pauli network's trials are no longer
reduced to their own cheapest before being compared with the direct
compilation. `through_pauli_network` returns every trial's result, in
order, and `paying_network` keeps the first of the cheapest of those
that need no more two-qubit gates and no more gates in all than the
direct candidates (`network_pays`, unchanged). `network_trials` adds four
trials (seeds 0–3) to the eight there were (seeds 0–7), synthesized with
`FRUGAL_LOCAL_PENALTY` (`core/src/pauli_network.rs`, 0.1): what the
greedy choice charges, in weights, for each single-qubit Clifford a
candidate CX needs. The usual `LOCAL_PENALTY` (0.01) only breaks ties.
`synthesize_seeded` takes the penalty as an argument.

New tests: `networks_are_judged_one_by_one` (the trial with the fewest
CXs that does not pay is passed over for one that does),
`both_penalties_are_tried_and_only_add_candidates`; the serial reference
of `candidates_side_by_side_choose_as_one_after_the_other` now runs every
trial one after the other and applies the new rule, and
`random_circuits_keep_their_unitary` checks both penalties.

**Why.** The log's standing item: on the 32-qubit Hamiltonian the
network saved 10% of the CXs but needed 26 gates too many to be kept.
Printing every trial showed two things. The trial kept for comparison was
the one with the fewest CXs (663 / 1581), so others that might pay were
never looked at; and none of the eight did (the fewest gates in all was
1581, against the direct 1555). A sweep of the local penalty (0.01, 0.05,
0.1, 0.2, 0.3, 0.5) on that circuit gave the best trials 663 / 1581, 666
/ 1563, 687 / 1511, 716 / 1500, 744 / 1472 and 823 / 1483: beyond 0.1
the CXs climb back to the direct compilation's. Both are changes to the
selection among candidates we already make, verifiable end to end today,
and a superset of the old trials, so no circuit can get worse. pytket's
GreedyPauliSimp has no such charge (its knobs are the discount rate,
depth weight, lookahead and trials), so this is ours. Portfolios were
compared on the benchmark's three Hamiltonians and 12 others (seeds
3000+, 8–32 qubits), then checked on 15 held out (seeds 5000+, 6–24
qubits): four trials at 0.01 and four at 0.1 (no more threads) gave 0.987×
the CXs on the tuning set but was worse on one tuning circuit and five
held-out ones; eight and eight was never worse, 0.979× on both sets, at
1.6–2.2× the time; eight and four was never worse, at 1.2–1.3×.

**Measured** (4-vCPU Xeon VM, complete runs of `run.py`, `device.py` and
`plugin.py`). All-to-all, two-qubit / all gates: Hamiltonian 8 107 / 293 →
101 / 286, 32 740 / 1555 → 687 / 1511 (Qiskit 740 / 1789, pytket 737 /
1681); every other row identical. Geometric means over Qiskit 0.938× →
0.933× the two-qubit gates, 0.877× → 0.875× all gates. On the chip, the
8-qubit Hamiltonian with CX 349 / 713 → 345 / 648, with CZ 331 / 1333 →
310 / 1119; no other row changed; QCC alone 0.935× / 0.825× → 0.934× /
0.820× with CX, 0.908× / 0.848× → 0.904× / 0.839× with CZ. Inside Qiskit
nothing changed (the plugin does not try networks). The 15 held-out
Hamiltonians: 0.981× the CXs and 0.980× the gates of 0.4.6, worse on
none. `check.py` passes. Time, 0.4.6 and 0.4.7 side by side, medians of
five, two rounds: Hamiltonian 8 29–30 → 40–44 ms, 16 130–133 → 173–180
ms, 32 325–335 → 402–451 ms (1.3–1.5×); other circuits try no network
and run the same code.

**Learned.** With the usual penalty, the Clifford left at the end of the
32-qubit network costs 29–47 CXs and 235–299 gates before optimization,
and the rotations 1472–1602 gates; the networks of the larger penalty
have 1617–1735 gates before optimization, against 1716–1848. Seeded noise and the
penalty find different networks: no portfolio of eight trials was
better everywhere.

**Next.**
1. Quantum Volume all-to-all speed (243 ms against Qiskit's 49 at 32
   qubits in this run): the serial passes after synthesis.
2. The Hamiltonian's compile time (now 0.4 s at 32 qubits against
   Qiskit's 15 ms): cheaper trials (stopping a trial's pipeline once it
   cannot pay), or a cheaper greedy step.
3. The few extra two-qubit gates against Qiskit on the chip (Quantum
   Volume 8: 111 against 108).

## 2026-11-06 (second change, asked for by the owner) — Shared syntheses and cheaper gathering, the same outputs

**What.** Two exact speed-ups of the pipeline that the CZ rounds repeat:

- `core/src/lower.rs`: `SynthesisCache`, the two-qubit syntheses of one
  compilation, shared by all its lowerings (one map per set of native
  gates, behind a mutex, on every thread). `Lowering::sharing` finds and
  leaves syntheses there; `prefetch_two_qubit` synthesizes only what no
  pipeline of the compilation has. `compile`, `compile_relabelling` and
  `compile_for_device` create one per call, and `passes.rs` and
  `device.rs` pass it down to every pipeline.
- `core/src/commuting_blocks.rs`: `gather_at` reorders indices into the
  circuit instead of copies of the region between the two blocks (it
  copied the region, every qubit's gates, once per cut and split tried),
  and remembers each range's CX count instead of working it out again for
  each cut; `crosses` reads the gates it crosses by reference.
  `pair_matrix` (`passes.rs`) is `Block::matrix` on borrowed gates, with
  the same products in the same order.

New test: `lower::tests::shared_syntheses_are_found_by_other_lowerings_
and_change_nothing`; `candidates_side_by_side_choose_as_one_after_the_
other` now runs every candidate through one cache shared with the earlier
circuits, against a reference that gives each pipeline its own.

**Why.** Picking up "Next" item 1 of the entry below (the CZ rounds' time).
Instrumenting the 16-qubit Hamiltonian on the CZ chip showed eight rounds
of 150–270 ms each, with the four cores already saturated, so only less
work helps. Two thirds of the syntheses (17 374 of 26 156, keyed by
unitary bits, objective and native gates) repeated one an earlier
pipeline or round had done, as most blocks come back unchanged, and
`gather_blocks` was the largest serial stage (17–145 ms per pipeline,
more than block resynthesis). Both can be made cheaper without changing
any choice. A one-pass CX → CZ translation would remove rounds, but it
changes outputs, and the measurements show why it is hard: the
through-CX candidate wins nearly every round on CZs while leaving 100–600
more single-qubit gates than the direct one (round 7: 1292 / 4572 against
1293 / 4451).

**Measured** (4-vCPU Xeon VM). Outputs: the 112 benchmark compilations
are byte-identical to 0.4.4's (and so to 0.4.5's); `check.py` passes;
`device.py`, a complete run, records the same gate counts, every state
checked. Time, 0.4.5 and 0.4.6 side by side, best of two medians of five,
geometric mean: on the chip with CZ 0.83× (0.65–0.99×), inside Qiskit
with CZ 0.82×, all-to-all with CZ 0.94×; with CX 0.96× on the chip, 1.00×
inside Qiskit and 1.02× all-to-all, noise (0.8–1.4× on circuits of a few
milliseconds). The 16-qubit Hamiltonian: on the CZ chip 2137 → 1601 ms,
inside Qiskit with CZ 1376 → 972 ms; Clifford+T 16 on the CZ chip 716 →
462 ms. `gather_blocks` there: 17–145 ms → 4–45 ms per pipeline.
`device.py` medians: QCC 41 ms with CX and 113 ms with CZ, Qiskit 39 and
38 ms in the same run.

**Learned.** With every core busy, the critical path does not set the
time: work removed anywhere counts. Sharing also lets later rounds skip
the threads of `prefetch_two_qubit` when fewer than 32 syntheses remain.

**Next.**
1. The rounds still cost most of the CZ time (1.6 s against Qiskit's
   0.15 s): the through-CX candidate's extra single-qubit gates (100–600
   per round) are what the next round folds; a better translation back to
   CZ, or a final direct pass kept only when cheaper, could end the rounds
   sooner without a worse result.
2. Quantum Volume all-to-all (142 ms against Qiskit's 61 at 32 qubits):
   the serial passes, where instructions allocate four times per copy.
3. Quantum Volume 8 on the chip: 111 two-qubit gates against Qiskit's 108.

## 2026-11-06 — Independent candidates side by side, the same outputs

**What.** `core/src/parallel.rs` (new): `map_ordered` and `join` run
independent work on all cores and return results in the order of the work,
on one thread on WebAssembly. Used in:

- `core/src/routing/sabre.rs`: `layout_and_route`'s trials (32 on the chip),
  each seeded by its number alone, side by side; the first of the cheapest
  is kept, as in the loop.
- `core/src/passes.rs`: `direct_candidates` compiles the direct pipeline and
  the one through CX side by side (for CZ-like devices, in every round of
  `until_fixed_point`); the Pauli networks are synthesized beside them and
  judged after (`beside_network`); below 64 instructions, one thread.
- `core/src/device.rs`: `with_and_without_network` gives the logical
  candidates to route from one compilation of the direct candidates (they
  were compiled twice, once in each `best_pipeline` call), and each
  candidate is routed and optimized on the chip on its own thread.

New tests: `parallel::tests` (order of results), `trials_side_by_side_
choose_as_one_after_the_other` (SABRE against a loop of single trials) and
`candidates_side_by_side_choose_as_one_after_the_other` (`best_pipeline`
and the device candidates against a serial reference of the old choice,
for CX and CZ, with and without networks).

**Why.** 0.4.4's repeated rounds made QCC 1.65× slower on the CZ chip (the
16-qubit Hamiltonian 4.4 s, against Qiskit's 0.13 s), and the log's
speed items were all open. A profile of the device path showed it serial
where the work is independent: 32 SABRE trials one after another, the
direct and through-CX pipelines one after another in each round, and the
direct candidates compiled twice. Running them side by side is exact (the
choice among candidates is unchanged by construction and checked byte for
byte), and could be verified end to end today; making the CZ translation
good in one pass, the other way to win the time back, changes outputs and
was tried before with mixed results (2026-11-03). On Quantum Volume 32
all-to-all, by contrast, the Rust core takes about 100 ms of QCC's 142 and
most of it is serial passes (resynthesis lookups, merging, phase folding),
so the Python wrapper (about 20 ms there, not a third) was left for later.

**Measured** (4-vCPU Xeon VM). Outputs: the 112 compilations of the
benchmark circuits (all-to-all, on the heavy-hex chip and inside Qiskit,
CX and CZ, 8–32 qubits) give byte-identical QASM and permutations against
0.4.4; `check.py` passes; `device.py`, a complete run, records the same
gate counts and every state checked (only its times changed). Time, 0.4.4
and 0.4.5 side by side, best of two medians of five, geometric mean: on
the chip 0.51× with CX (0.29–1.08×) and 0.62× with CZ (0.39–1.05×); inside
Qiskit with CZ 0.82×; all-to-all with CZ 0.85×. The 16-qubit Hamiltonian
on the chip: 1377 → 566 ms with CX, 5583 → 2205 ms with CZ; Clifford+T 16
with CZ 1291 → 700 ms; Quantum Volume 16 with CX 463 → 180 ms. Where
nothing new runs side by side (CX all-to-all, CX inside Qiskit): 1.04× and
1.00× in this run, 0.98× and 0.98× in the one before, individual circuits
of a few milliseconds 0.8–1.4×: noise. `device.py`'s median for QCC with
CX went from 64 ms (0.4.4's run) to 38 ms, Qiskit's 36 ms in the same run.

**Learned.** A thread costs a few tenths of a millisecond: the first
version spawned one for every compilation and made 8-qubit GHZ with CZ
1.2× slower (1.5 → 1.7 ms); hence the 64-instruction threshold and no
thread when there is no second candidate. Under valgrind threads run one
at a time, so callgrind's per-thread totals say nothing about balance.

**Next.**
1. The CZ rounds themselves (the 16-qubit Hamiltonian on the chip is still
   2.2 s against Qiskit's 0.15 s): a one-pass CX → CZ translation, or
   rounds that resynthesize only the blocks the last round changed.
2. Quantum Volume all-to-all (QCC 142 ms against Qiskit's 61 ms at 32
   qubits): the serial passes after synthesis — instructions own a `String`
   and three `Vec`s, so every copy and relabel allocates four times (malloc
   and free are a fifth of the instructions); the two measures' block
   resyntheses and merged costs could run side by side.
3. The few extra two-qubit gates against Qiskit on the chip (Quantum
   Volume 8: 111 against 108).

## 2026-11-03 — Optimizing again until the circuit stops shrinking, for CZ devices

**What.** `core/src/passes.rs`: `until_fixed_point` runs the default
pipeline again on its own result, without relabelling or new pairs, while
that strictly lowers (two-qubit gates, all gates), at most eight times. It
runs only for devices whose two-qubit gate is CX-like but not CX (CZ, ECR,
…), in `compile` (the Qiskit plugin), `compile_relabelling` (all-to-all) and
on the physical qubits after routing (`on_the_same_pairs`, used by
`core/src/device.rs`). New test:
`repeated_rounds_keep_the_unitary_and_the_pairs_and_never_cost_more`
(CZ and ECR circuits get cheaper on some inputs and are never worse, keep
the unitary and the input's pairs; for CX the result is identical).

**Why.** On the heavy-hex chip with CZ, Qiskit routing plus QCC's
optimization left fewer gates in all than QCC alone (16-qubit Hamiltonian
4509 against 5372; Clifford+T 2579 against 2880, where QCC alone used more
gates than Qiskit's 2844). Running QCC's pipeline a second time on its own
CZ output removed more gates every time (5372 → 5205 → … → 4633 over five
rounds), while on CX it was already a fixed point. The cause: for CZ the
pipeline reads the circuit through CX (each CZ as H·CX·H, the H on its
second qubit) for phase folding, then translates it back, and every
translation leaves the single-qubit gates where the next reading finds
parities and blocks the last did not. Repeating to a fixed point is what
Qiskit's level 3 does with its optimization loop (`DoWhileController` with
`MinimumPoint` on depth and size), and alternating cheap rewrites with
resynthesis until nothing improves is the core of GUOQ (Xu, Molavi, Tannu
and Albarghouthi, "Optimizing Quantum Circuits, Fast and Slow", ASPLOS 2025,
arXiv:2411.04104). It is exact (every round is the verified pipeline, kept
only when strictly cheaper) and could be verified end to end today. A
per-CZ choice of which side gets the H (the side whose neighbouring runs
become diagonal, or need fewest √X) was tried as a Python prototype: it won
on some circuits (Hamiltonian 8 on the chip, 1352 → 1206 gates after a
second pass) and lost on others, so it was left out.

**Measured** (4-vCPU Xeon VM; `device.py` and `plugin.py` complete runs,
every output's state checked; `check.py` passes; `run.py` is CX only and
unaffected, its outputs identical by construction and by `check.py`).
With CZ every row got cheaper or stayed the same, none worse; with CX none
changed. On the chip, QCC alone with CZ: 0.917× → 0.908× Qiskit's two-qubit
gates, 0.868× → 0.848× its gates; Hamiltonian 16 1379 / 5372 → 1292 / 4572
(Qiskit 1477 / 5357), Hamiltonian 8 336 / 1378 → 331 / 1333, Clifford+T 16
722 / 2880 → 708 / 2790, 8 234 / 865 → 227 / 826, QAOA 16 182 / 726 →
180 / 704. Inside Qiskit with CZ: 0.947× / 0.861× → 0.936× / 0.837×;
QFT 16 1312 → 1244 gates, Hamiltonian 16 1405 / 4509 → 1312 / 4416,
Heisenberg 8 430 → 391. Time, against 0.4.3 side by side (median of five):
CZ on the chip 1.65× in geometric mean (Hamiltonian 16 1.69 s → 4.37 s,
Clifford+T 16 0.45 → 0.97 s), the plugin's pipeline with CZ 2.5× (tens of
milliseconds; Heisenberg 16 18 → 108 ms); CX 0.99× / 1.00×.

**Learned.** The first physical pass on the CZ Hamiltonian chose the
through-CX candidate for 37 fewer CZs at the price of 768 more gates (1379 /
5372 against the direct 1416 / 4604); the rounds then descend slowly, and
the eighth still removed CZs. The translation from CX back to CZ is where
the single-qubit gates come from.

**Next.**
1. Make the CX → CZ translation good enough in one pass, so the rounds are
   not needed (and the time comes back): choose each CZ's H side with the
   whole run in view, or resynthesize the runs around each CZ with the
   neighbouring blocks.
2. The Python wrapper: about a third of Quantum Volume 32's time (see
   2026-11-01).
3. Win back the adder's gathered CNOTs without the single-qubit cost (see
   2026-10-30).

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
