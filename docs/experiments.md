# Experiments: what worked, what didn't, what's next

This file is the compiler's memory across improvement runs. Before choosing
what to work on, read **Next** and search **Didn't work** and **Considered**
for the idea, its file and its source. Don't retry something listed there
unless you say what is different this time. At the end of every run, kept or
not, update all four sections in the same commit as the change; if nothing
ships, commit this file alone.

`docs/improvement-log.md` has the full story of each shipped change; this
file is the short index, including everything that was tried and dropped.
Gate counts are two-qubit / all gates. "Chip" is the 19-qubit heavy-hex
benchmark device. Times are from the run's own machine (a 4-vCPU VM unless
stated) and only comparable within that run.

## Next

In priority order. Replace this list at the end of each run.

1. **CZ device time, by doing less per round.** The 16-qubit Hamiltonian
   on the chip still takes 2.57 s (0.4.12; Qiskit 0.11–0.15 s). Pruning
   candidates does not work (2026-11-17, **Didn't work**), and the
   profile of one candidate's rounds is spread: synthesis misses 28%
   (`dress`, `kak`, `cx_count`), gathering 11% (5.3% of it `cx_count`'s
   KAKs), M4 products 10%, rotation merging 5.5% after 0.4.12. The big
   lever left is incremental rounds: rerun only the parts of the circuit
   the last round changed (each round changes a few blocks; the cache
   already shows it, 66% of syntheses repeat). Needs a way to bound how
   far a change can reach through phase folding and gathering.
2. **Carrying, beyond the strict rule** (2026-11-15). Applied to the
   final output and judged by CX alone, carrying would take chip CX
   Clifford+T-16 from 733 / 1209 to 704 / 1251 and QAOA-16 all-to-all from
   96 / 317 to 90 / 340; the no-more-of-either rule refuses both (and
   `check.py` would fail on the all-to-all one). A carry that leaves fewer
   single-qubit gates (choosing the split U = U₂·U₁ by the merged blocks'
   synthesized gate counts, not only cut at gate boundaries; or carrying
   before `shift_through_blocks` and letting it absorb the R gates) could
   win part of it within the rule. Per-move judging by synthesized gate
   counts (below, **Didn't work**) was not enough.
3. **More routed candidates for quality** (2026-11-17 data). The first
   pass ranks the routed candidates badly (on one held-out circuit the final
   winner started more than 7% behind), so optimizing more routings than the
   current 2–4 could pay in gates; it was small at 3.5× the time in
   2026-11-09 ("Top-k routings"). Retry only once rounds are much
   cheaper (item 1).
4. **Quantum Volume speed, all-to-all.** QV-32 96 ms in the 2026-11-17
   run, of which the core 78 ms (the wrapper is now a quarter). Synthesis
   (`TwoQubitSynth::dress`, `Shape::of`, `M4` products, `hypot`).
5. **Hamiltonian compile time, continued.** 162 ms at 32 qubits against
   Qiskit's 18 in the 2026-11-15 run. Left: 13 full pipelines (one per
   network trial, plus the direct one) and the synthesis's `layers()`.
   Stopping a trial's pipeline once it cannot beat the direct compilation
   is not exact; check on held-out Hamiltonians.
6. **Routing quality beyond more trials.** QV-8 on the chip still uses 111
   CNOTs (Qiskit 108 or 111 depending on its seed). Untried from
   LightSABRE (Zou et al., arXiv:2409.08368): its relative scoring;
   routing whole two-qubit blocks as one node.
7. **The adder's gathered CNOTs** (37/85/181 at 8/16/32 qubits, against
   40/92/196 kept). They cost 108/252/540 gates against 105/245/525, so the
   keep rule declines them. Win them back by choosing the outer gates of the
   SWAP-elided blocks with their neighbours in view.
8. **Pauli-network judging after carrying** (2026-11-15, not observed).
   `network_pays` compares a network with the direct candidates after
   both were carried; carrying could in principle make the direct one
   cheaper in all gates while the network keeps fewer CXs, flipping the
   choice to more CXs. None of 30 Hamiltonian compilations showed it;
   fix the same way as gathering (judge on costs without carrying) if one
   does.

## Worked

Newest first. Details are in `improvement-log.md` under the same date.

| Date | Version | Change | Effect |
|---|---|---|---|
| 2026-11-17 | 0.4.12 | One update of the Clifford map per gate in rotation merging and the network's extraction (`apply_all_to_map`, `pauli_network.rs`: √X, X, Y, CZ, SWAP and π/2-multiple Z rotations were 2–6 updates each); `read` returns its steps inline instead of a `Vec` | 220/220 outputs byte-identical (incl. 81 held-out). One CZ chip candidate of Hamiltonian-16: 15.30 → 14.74 G instructions, rotation merging 1.37 → 0.81 G. Paired: all-to-all 1.00× CX, 0.974× CZ; chip 1.001× CX, 0.983× CZ; CZ Hamiltonian-16 2.77 → 2.57 s |
| 2026-11-15 | 0.4.11 | Carrying gates between blocks on the same pair (`carry_through`, `commuting_blocks.rs`): the end of a block (from one of its last 3 two-qubit gates) forward, or the start of one back, through single-qubit gates by conjugation (g D g†) and through blockers that commute exactly, where the two blocks then need fewer CX-like gates; only on final results, after resynthesis, kept only if no more of either; gathering judged on costs without carrying; CZ rounds plain first, then carrying (`Last` in `passes.rs`) | All-to-all Clifford+T 16 150/505 → 148/499 (pytket 148/495), 32 310/1034 → 308/1030; geometric mean 0.933×/0.875× → 0.932×/0.874×. Chip CX: QAOA-8 94 → 92, adder-16 148/405 → 146/402, Hamiltonian-8 345/648 → 335/639, -16 1336/2100 → 1321/2099; CZ: Clifford+T-16 692 → 691, Hamiltonian-8 310/1119 → 307/1103, -16 1256/4454 → 1231/4208; inside Qiskit 5 rows fewer 2q. Held-out 82: fewer 2q on 22, more on none. Time 1.05× all-to-all, 1.07× chip (CZ Hamiltonian-16 1.46 → 2.44 s) |
| 2026-11-13 | 0.4.10 | Pauli-network scoring from per-candidate tables (`candidate_tables`), nine scores summed in one pass in the old order, weights read once per step, √X conjugated in one step (`pauli_network.rs`); per-pass KAK memo and skipping free sides of lone CXs (`symmetric_blocks.rs`) | 142/142 outputs byte-identical (incl. 30 held-out Hamiltonian compilations). Ham-32 core 14.4 → 3.4 G instructions (synthesis 12.0 → 1.4 G). All-to-all Hamiltonian 8/16/32 0.59×/0.37×/0.32× the time with CX (340 → 110 ms at 32), 0.71×/0.52×/0.42× with CZ; others 0.94×; chip 0.96× (CZ Ham-16 1.45 → 1.37 s) |
| 2026-11-11 | 0.4.9 | Instructions own no heap memory (`ir.rs`: `Name` borrows known gate names, `Inline<T, 3>` qubits and parameters); FxHash for internal maps (`hashing.rs`); Python's cyclic GC paused while instructions cross to and from the core (`paused_gc`, `native.py`, also the Qiskit plugin) | 160/160 outputs byte-identical. All-to-all 0.85× the time (geometric mean, side by side); QV-32 149 → 92 ms (core alone 93 → 59 ms, 1.09 → 0.91 G instructions, allocation 27% → 10%); chip 0.93× CX, 0.83× CZ (CZ Hamiltonian-16 1.79 → 1.42 s). Faster than Qiskit on 20 of 24 circuits (was 18) |
| 2026-11-09 | 0.4.8 | Faster SABRE with the same choices (`routing/sabre.rs`: shared `Dag`, scratch-space extended set, `Trace` written out only if kept, `Savings` memo, flat distances); device routing 16 swap trials per layout instead of 4, and the 4-trial routing optimized beside it (`layout_and_route_and_first`, `device.rs`) | Routing 2.1–2.9× fewer instructions, same routings (112/112 outputs byte-identical with old options). Chip: QV-16 822 → 783 CX/CZ, Clifford+T-16 752 → 733 CX, 708 → 692 CZ, Hamiltonian-16 1368 → 1336 CX, 1292 → 1256 CZ, QAOA-16 CX 180 → 179; no row worse. Held-out 40: fewer 2q gates on 19, more on none. Time 0.86× CX, 0.97× CZ (CZ Hamiltonian-16 1.37×) |
| 2026-11-07 | 0.4.7 | Pauli-network trials judged one by one against the direct compilation; four extra trials at local-Clifford penalty 0.1 (`FRUGAL_LOCAL_PENALTY`) beside the eight at 0.01 | Hamiltonian-32 740/1555 → 687/1511; Hamiltonian-8 107/293 → 101/286; chip CZ Hamiltonian-8 331/1333 → 310/1119; geometric mean vs Qiskit 0.938× → 0.933× two-qubit gates; 15 held-out Hamiltonians 0.981× CNOTs, worse on none; Hamiltonians compile 1.3–1.5× slower |
| 2026-11-06 | 0.4.6 | Compile-wide cache of two-qubit syntheses per native set (`SynthesisCache`, 66% of requests repeat); `gather_blocks` reorders indices instead of copying | Same outputs. Chip CZ 0.83× the time; chip CZ Hamiltonian-16 2137 → 1601 ms |
| 2026-11-06 | 0.4.5 | Independent candidates side by side (`parallel.rs`): SABRE's 32 layout trials, direct and through-CX pipelines, Pauli networks, device candidates; 64-instruction threshold before spawning | 112/112 outputs byte-identical. Chip 0.51× the time with CX, 0.62× with CZ |
| 2026-11-03 | 0.4.4 | `until_fixed_point`: for CX-like non-CX devices (CZ, ECR), rerun the pipeline while it strictly shrinks, at most 8 rounds | Chip CZ two-qubit 0.917× → 0.908× Qiskit, all gates 0.868× → 0.848×; Hamiltonian-16 chip CZ 1379/5372 → 1292/4572. CZ compiles 1.65× slower on the chip, 2.5× in the plugin |
| 2026-11-01 | 0.4.3 | Exact speed-ups: `count_1q` without Euler angles, sparse parity lists in phase folding, rotation merging only over a generator's support, one `Plan` per block in resynthesis | Byte-identical outputs; 0.73× the time in geometric mean; QV-32 230 → 155 ms |
| 2026-10-30 | 0.4.2 | Gather same-pair two-qubit blocks through commuting gates (`commuting_blocks.rs`, Fagan & Duncan arXiv:1901.10114); single-qubit trades through lone CX blocks in the last pipeline only | Clifford+T 154/525 → 150/505 at 16 qubits; chip CZ 0.935×/0.878× → 0.917×/0.868×; adder about 2.5× slower to compile |
| 2026-10-29 | 0.4.1 | Pauli networks with lazy local frames, CX orientation by local cost, √X basis changes; eight seeded trials (noise 0.05 × seed) | Hamiltonian-16 325/1285 → 313/792 (Qiskit 393/972, pytket 390/887); geometric mean 0.940×/0.882× |
| 2026-10-19 | 0.4.0 | Pauli networks (Rustiq / GreedyPauliSimp style), rotation merging through Cliffords, Hamiltonian benchmark | 0.942×/0.916× Qiskit O3 |
| 2026-10-19 | 0.2.2–0.3.0 | Phase folding, in-context two-qubit synthesis, SWAP elision, symmetric blocks, CZ via CX, layout and SABRE routing with NASSC costs | See `improvement-log.md` and git history |

## Didn't work

Don't retry these as they were tried. Each line says why. "Retry if" says
what new fact would justify another attempt.

**Pauli networks (`pauli_network.rs`, `passes.rs`)**
- **Weighted keep rule (a CX counts as 10 gates)** (2026-10-29). Failed Pareto selection at 16 qubits and regressed total gates at 32. Replaced by "no more two-qubit gates and no more gates in all".
- **Tracking each wire's pending local Clifford to score merges** (2026-10-29). Didn't beat counting local gates.
- **Seed noise ≥ 0.8** (2026-10-29). Loses to the deterministic run.
- **Environment-variable toggles in the hot loop** (2026-10-29). They slowed synthesis about 2×. Use them only in throwaway prototypes and remove them before timing.
- **Local-Clifford penalty sweep** (2026-11-07, 32-qubit Hamiltonian, direct 740/1555). Best trial per penalty: 0.01 → 663/1581, 0.05 → 666/1563, 0.1 → 687/1511 (kept), 0.2 → 716/1500, 0.3 → 744/1472, 0.5 → 823/1483. Above 0.1 the CNOTs climb back.
- **Trial portfolios** (2026-11-07, tuning set of 15 plus 15 held out).
  - 4 at 0.01 + 4 at 0.1 was worse on 6 circuits.
  - 8 + 8 was never worse but 1.6–2.2× slower.
  - 8 + 4 shipped: never worse, 1.2–1.3× the time.
  - No portfolio of eight was better everywhere: seeded noise and the penalty find different networks.
- **Pauli network on Clifford+T circuits** (2026-10-30). Much worse: direct vs network, CNOTs/gates: 8 qubits 77/287 vs 103/297; 16 qubits 155/521 vs 255/681; 32 qubits 316/1065 vs 906/2217. The network is only tried on circuits with three or more CXs per rotation; keep it that way.

**Resynthesis of larger regions (prototypes on QCC's output, 2026-11-13)**
- **Phase-polynomial resynthesis of CNOT+Rz regions** (GraySynth, Amy, Azimzadeh and Mosca, arXiv:1712.01859, then PMH for the residual linear map; Python prototype on QCC's all-to-all outputs, regions grown greedily through cx/rz with single-qubit gates allowed in front). Saved 0 CNOTs on all 21 circuits; worse on most regions (Clifford+T-32: worse in 21 of 21, Hamiltonian-32: 70 of 97). Regions are small in QCC's output (√X everywhere; at most 4–28 CXs, QAOA's 96 aside, where it was also no better). Retry only on input circuits before phase folding, or with an optimal small-region CNOT synthesis.
- **pytket `ThreeQubitSquash` on QCC's output** (Next item 6 of 2026-11-11). All-to-all, without swaps: Clifford+T 16/32 150 → 148, 310 → 308, Hamiltonian-16 313 → 310, everything else unchanged; with swaps the adder takes the known gathered-CNOT trade (37/85/181) and Clifford+T-32 306. On the chip it looked large (Hamiltonian-16 1336 → 1233, -8 345 → 314, Clifford+T-16 733 → 704, QAOA-16 179 → 173) but put 323, 71, 3 and 2 CXs on uncoupled pairs: invalid. A coupling-aware three-qubit synthesis would be needed there. The all-to-all wins came from small structures (see Next item 4), not from blocks above 20 CXs.

**Gathering and single-qubit trades (`commuting_blocks.rs`, `symmetric_blocks.rs`)**
- **Moving a blocker out of the whole surrounding runs** (2026-10-30). The adder only went 40 → 39, because a √X earlier in a run stops the blocker. Cutting inside the runs and carrying third-qubit gates works.
- **Gathering without the keep rule** (2026-10-30). The adder reached 37/85/181 CNOTs but more gates in all. The pipeline now runs with and without gathering and keeps the cheaper.
- **CX trades before routing, or in the CX stage of CZ compilation** (2026-10-30). Chip CZ Clifford+T-8 cost 100 more gates, because phase folding lost its Z rotations. Trades run only in the last pipeline (`last`).
- **A guard keeping diagonal runs diagonal** (2026-10-30). No gain (CZ Clifford+T-8 still 236/965). Running the trades last was the fix.
- **Accepted weaknesses after 0.4.2** (2026-10-30). Some chip rows got slightly worse in one count from routing variance: CX Clifford+T-16, CZ adder-16, CZ Hamiltonian-8, CZ QFT-16. Remedy, not yet tried: route the logical circuit with and without gathering and keep the cheaper.

**Carrying gates between blocks (`commuting_blocks.rs` `carry_through`, `passes.rs`, 2026-11-15)**
- **Carrying with gathering, before resynthesis** ("placement A": `carry_through` right after `gather_blocks`, every pipeline). Chip CX Hamiltonian-16 1336 → 1373 CNOTs, CZ Clifford+T-8 227/826 → 229/838; total two-qubit gates over the 80 benchmark compilations unchanged. It changes the circuit before routing and the CZ rounds' trajectory. Carrying runs only on final results.
- **Carrying in every pipeline after resynthesis, including before routing** ("B" without `last`): CZ rows worse (Clifford+T-16 692 → 695, Hamiltonian-16 1256 → 1263 in one variant), chip CX Hamiltonian-16 sometimes worse. Restricting it to `last` pipelines fixed the CX rows.
- **Keeping the carried circuit when lexicographically cheaper** (fewer CXs first, "L"). The largest wins (chip CX Clifford+T-16 733/1209 → 704/1251 with all moves; 716/1228 in one variant; all-to-all Hamiltonian-16 313/790 → 310/801) come with more gates in all, which the pipeline's no-more-of-either rule and `check.py` (all-to-all) refuse. Retry only with a decision on that rule (see Next 2).
- **Per-move judging by synthesized gate counts** (each move also had to leave the two blocks, synthesized with `Objective::InCircuit`, with no more gates in all, or no more inner single-qubit gates). Not better than judging by CX alone under the global rule: all gates 1321/2099 vs 1324/2094 and 1328/2086 (inner) on chip CX Hamiltonian-16, but fewer CNOTs saved over the benchmark (−40 and −58 against −67), and a synthesis per candidate. Dropped for the simpler CX-only test.
- **Carrying in the CZ rounds from the start.** CZ Clifford+T-16 on the chip 692/2677 → 695/2688: a different trajectory of rounds ends at a worse fixed point. Shipped: rounds without carrying until they stop, then rounds with carrying (never worse).
- **More blockers (8 instead of 3)**: CZ QAOA-16 180/704 → 175/713, Hamiltonian-16 1231/4208 → 1226/4236, 1.19× the time; mixed. **One carried two-qubit gate instead of up to 3**: QAOA-8 chip 92 → 94. **Up to 6**: no change. Kept 3 and 3.
- **Ruled out while debugging**: an all-to-all Clifford+T-16 result of 152/504 (more CXs than 0.4.10's 150/505) came from the gathering keep rule, not carrying itself: carrying made the circuit without gathering cheaper in all gates, so the gathered one no longer dominated it and the rule fell back to the one with more CXs. Fixed by judging gathering on the costs without carrying.
- **Skipping the carrying rounds when the result has nothing to carry** (shipped as a speed-up) lost one row's 7 gates (all-to-all CZ Heisenberg-32 372/1789 back to 0.4.10's 372/1796).

**Pruning routed candidates on a chip (`device.rs`, 2026-11-17, Next item 1 of 2026-11-15)**
- **Rounds only for candidates whose first pass is within a margin of the best** (2-qubit count after `first_pass`, before `until_fixed_point`). Traced 97 CZ compilations (the benchmark's 16 and 27 held-out sizes on the chip, 31 other-seed circuits on the chip, 23 on a 4×5 grid), 61 with 2–4 candidates, 1518 rounds. Results changed by margin: 0% 19, 1% 13, 2% 6, 3% 3, 4% 3, 5% 1, 7% 1, 10% 0; at 10% only 145 of 1518 rounds (10%) would go. Example: held-out Hamiltonian-14 on the chip, the winner (965 CZ) started 2.1% behind (1060 vs 1038) and the other ended at 992. The first pass ranks candidates poorly because the rounds take 2–11% off, unevenly. Retry only with a predictor better than the cost so far (none found).
- **Lockstep pruning after every round** (same data, same margin on the current cost): 3% 4 changed, 5% 1, 10% 0 with 15% of rounds saved. No better.
- **Dropping one of the two pipelines in the rounds** (direct, or through CX). Over the 43 chip CZ compilations the direct one was cheapest in 368 rounds (365 strictly), the one through CX in 139: both are needed.
- **A memo of `cx_count` by the unitary's bits** (gathering's `cx_cost`, 5.3% of one CZ candidate's instructions). 14 966 hits against 19 159 misses on the CZ Hamiltonian-16 candidate: at most 2.3% saved, not worth a cache threaded through `gather_blocks` or a global one.

**CZ devices (`passes.rs`, `until_fixed_point`)**
- **Rounds without the last-pass single-qubit shift, then one final full pass ("mode 2")** (2026-11-03). Hamiltonian-16 4430 gates in 4896 ms against mode 1's 4572 in 4445 ms. But it was worse on Clifford+T-8 (229 vs 227), QAOA-16 (713 vs 704), QAOA-8 (363 vs 360), Hamiltonian-8 (1334 vs 1333) and slightly worse on CX rows. Mode 1 shipped.
- **Per-CZ H-side orientation** (2026-11-03, Python prototype on device output plus one recompile). Choose which qubit of each CZ takes the H, scored by "neighbouring runs become diagonal" or by √X count. Gates against a plain second pass:
  - Hamiltonian-8: 1352 → 1219 / 1206 (330 CZ against 335).
  - QFT-16: 1421 → 1357 / 1412.
  - Hamiltonian-16: 5205 → 5195 / 5185.
  - Clifford+T-16: 2813 → 2815 / 2811.
  - QAOA-16: 708 → 705.
  - Adder-16: unchanged.
  - The all-to-all variant was worse: QAOA-16 452 vs 440.
  - Left out as mixed. Retry if implemented in Rust and combined with the rounds, choosing the side with the whole run in view; it was never tried that way.
- **CX repeated rounds** (2026-11-03). Rerunning `compile` on CX output still shaves 1–4 gates (Clifford+T-16 500 → 496), never two-qubit gates. Not worth the time, so the rounds skip CX.

**Routing (`routing/sabre.rs`, `device.rs`)**
- **More layout trials** (2026-11-09, 64 or 128 instead of 32, 4 swap trials). QV-8 stayed 111; 128 helped QV-16 (813) and the adder (144) but at 3× the routing time, and 64 made the CZ Hamiltonian-8 worse (332 vs 310). More swap trials per layout did more per unit of time.
- **Swap trials 8 or 32 alone** (2026-11-09). 8: CX Hamiltonian-16 1368 → 1382, Clifford+T-16 752 → 756. 32: CZ Hamiltonian-16 4572 → 4913 gates. 16 alone was never worse on the 32 chip rows but was on 3 of 40 held-out circuits (Clifford+T-10 CX 342 → 350, Hamiltonian-14 CX 1068 → 1082, Hamiltonian-10 CZ 465 → 483): the SWAPs' CNOT estimate is not the optimized count. Shipped: 16, with the 4-trial routing optimized too.
- **Top-k routings by estimate, each optimized** (2026-11-09, 4 swap trials, k = 4 or 8). Small wins (QFT-8 CX 80 → 77, adder-16 148 → 146, CZ QAOA-8 93 → 89), QV-8 unchanged, CZ time up to 3.5× (Hamiltonian-16 1.9 → 6.6 s at k = 8). Retry if the physical pipeline gets much cheaper.
- **Block-aware extended set** (2026-11-09: a two-qubit gate continuing a block on the same pair not counted again in the lookahead, as if routing blocks). QV-16 822 → 807, Hamiltonian-16 1368 → 1348, QAOA-16 180 → 173, but QFT-16 353 → 438 and Clifford+T-16 752 → 762. Retry as part of a real block-level router with its own tuned weights.
- **merge_weight 0 or 1** (2026-11-09). 0 (plain SABRE): QFT-16 353 → 376, adder-16 148 → 170. 1: same as 0.5 on the rows tried. Kept 0.5.
- **Ruled out while debugging QV-8's 111**: the SWAP-cost estimate is exact there (72 logical CX + 39 = 111 for every routing printed), so the gap is the search, not the optimization after it; and Qiskit itself gives 108 or 111 with seeds 1–10.

**Speed (`lower.rs`, `parallel.rs`, `ir.rs`, `pauli_network.rs`, Python wrapper)**
- **Ruled out while profiling the Hamiltonian** (2026-11-13). The Next list's guesses (stop a trial early, cache the templates' KAK per Weyl point) were not where the time went: 83% was the network's candidate scoring. Summing the nine candidates' scores per Pauli in another order (grouped by letter pair) would be faster still but changes float rounding, so the order was kept.
- **Ruled out while profiling the wrapper** (2026-11-11). Building the raw tuples (5.5 ms) and the output `Instruction` dataclasses (5 ms) on QV-32 are not the wrapper's cost; the garbage collector was (15–20 ms, see **Worked** 0.4.9). The 2026-11-01/03/06 entries that blamed the dataclasses were misread.
- **A thread for every compilation** (2026-11-06). 8-qubit GHZ with CZ got 1.2× slower (1.5 → 1.7 ms): a thread costs a few tenths of a millisecond. There's now a 64-instruction threshold, and no thread without a second candidate.
- **Reading a profile as "no single hotspot"** (2026-10-30). The QV-32 profile looked spread out (resynthesis 44%). A callgrind look on 2026-11-01 found the serial passes were the bottleneck and the block syntheses on 4 cores were not. Profile the main thread, not totals, and note that under valgrind threads run one at a time.

## Considered, not tried

- **Incremental CZ rounds** (2026-11-17). Rerun the pipeline only where the last round changed the circuit. Likely the largest lever on CZ chip time (see **Next** 1), but phase folding and gathering reach across the whole circuit, so the region a change can affect must be bounded first; not started.
- **External research skipped** (2026-11-17). This run went straight to Next item 1 (a measured weakness) and spent its time on the pruning data and the profile; no new paper was read.
- **Python wrapper on Quantum Volume** (2026-11-01, 2026-11-03, 2026-11-06). Done on 2026-11-11 in another way: most of it was the garbage collector, now paused (`paused_gc`). Left: `PassManager.run`'s circuit copy (3.8 ms on QV-32), building the tuples and dataclasses (about 5 ms each way), and pyo3's own conversion (**Next** item 2).
- **A global allocator (mimalloc, jemalloc) in the Python extension** (2026-11-11). Would cut the remaining 10% allocation cost without code changes, but adds a C dependency to every wheel target of the Release workflow and cannot serve the WebAssembly build; removing the allocations was exact and portable. Retry if allocation is again a large share after the remaining sources (phase folding's groups, `dress`, `candidate_angles`) are gone.
- **`smallvec` / `rustc-hash` crates** (2026-11-11). The core has a single dependency (`num-complex`); the two small types it needed (`Inline`, `FxHasher`) are about 250 lines with tests in `ir.rs` and `hashing.rs`.
- **A new synthesis method instead of exact speed-ups** (2026-11-01). Speed-ups were chosen because they can be verified byte for byte.
- **LightSABRE (Zou, Treinish, Hartman, Ivrii and Lishman, arXiv:2409.08368)** (2026-11-06, 2026-11-09). Read for 0.4.8: QCC already has its release valve and uses SABRE's nondeterminism (trials); 0.4.8 took its lesson that more trials pay once routing is cheap. Its relative scoring would change float ties, so the outputs could not be checked byte for byte; left for a run that aims at routing quality (Next item 4).
- **CZ-native single-qubit reduction** (2026-11-06). Searched; no decision recorded.
- **Splitting a block to move part of it past a gate** (2026-11-13; built 2026-11-15 as `carry_through`, splits at gate boundaries only). From bisecting where `ThreeQubitSquash` saves CNOTs on Clifford+T-16: CX(15,1) H(1) CX(1,7) CX(15,1) = CX(15,1)·CZ(15,1) after CX(1,7), H(1) first, as H·CX(15,1) = CZ(15,1)·H and CZ commutes with CX(1,7)'s control. Needs a search over splits U = U₂·U₁ (U₂ diagonal on the shared qubit in the middle gate's basis); not built this run because the gain measured all-to-all was 2–3 CNOTs and the chip gain is unknown without a prototype.
- **A coupling-aware three-qubit synthesis** (2026-11-13). The only way pytket-style three-qubit squashing could be used on a chip (heavy-hex triples are paths a–b–c). Large (a QSD-style decomposition restricted to two edges); not tried.
- **pytket `ThreeQubitSquash`** (2026-10-30). Measured 2026-11-13: see **Didn't work**, "Resynthesis of larger regions".
- **A continuous search for the split U = U₂·U₁** (2026-11-15). U₂ must commute with the gates between (for a CX control on the shared qubit: block-diagonal on it, a controlled gate with one CX at most), a 7-parameter family; minimizing the two blocks' CX counts over it would find splits that are not at gate boundaries. Not tried: the gate-boundary splits found the benchmark's cases, and the CX count is a step function, hard to optimize numerically (the Makhlin-invariant fit of `numeric.rs` would be the starting point).
- **ZX-calculus simplification and extraction of whole circuits** (PyZX's `full_reduce`; Duncan, Kissinger, Perdrix and van de Wetering, arXiv:1902.03178) (2026-11-15). Aimed at T-count and gate count on Clifford+T circuits; extracting a circuit from the reduced graph is not coupling-aware and can raise the CNOT count, so it would not serve the chip. Not measured on QCC's output; a quick upper-bound check like the 2026-11-13 `ThreeQubitSquash` one would tell.

## Method lessons

- **Timing on the VM.**
  - Compare builds only side by side, interleaved, in the same run, with medians of 9–11.
  - Short runs (a few ms) need hundreds of calls.
  - With every core busy, any removed work counts, not just the critical path.
- **Tests at tolerance edges.** A count test exactly at the tolerance (θ − π/2 = 1e-10) can go either way under rounding in both methods. The 2026-11-01 run changed the test's nudges to either side of the edge. Never weaken a test to pass; if an edge is genuinely ambiguous, say so in the log.
- **Outputs.** For exact speed-ups, compare all 112 benchmark compilations byte for byte before and after.
- **Python tests need `qevedo-simulator`** (2026-11-11). Without it (`pip install ".[dev]"` brings it), `test_native.py` and `test_synthesis.py` skip as whole modules and pytest reports 86 passed, 2 skipped instead of 219 passed.
- **Timing Python glue** (2026-11-11). Time it with `gc.disable()` as well as without: bulk object creation sets off collections that walk the whole heap (Qiskit included), and that, not the objects themselves, was most of the wrapper's cost. cProfile hides it inside whichever function happened to allocate.
- **Byte-for-byte outputs, quickly** (2026-11-11). A throwaway script compiling all 24 circuits all-to-all and 16 on the chip, with CX and CZ, QCC alone and through the plugin (160 outputs, about 25 s), written to one file before and after and compared with `cmp`; two runs of the old build gave the same file, so the outputs are deterministic. A second venv with the old wheel (`maturin build` from a `git archive` copy) gives paired end-to-end timings.
- **Profiling a Rust pass without Python** (2026-11-09). Dump the pass's input from the device path (a throwaway env var) and time it in a throwaway `core/examples/` binary under `valgrind --tool=callgrind`; instruction counts are stable where wall time on this VM varied ±20%. With `CARGO_PROFILE_RELEASE_DEBUG=true`, `callgrind_annotate --auto=yes` gives line-level counts.
- **Python string edits that silently miss** (2026-11-09). A `str.replace` on code that an earlier experiment had changed did nothing, and the default stayed at 4 swap trials until a test caught it. Assert that the old text is there before replacing.
- **Tests that must not be vacuous.** A "never worse" test should also assert that the two sides differed somewhere (`the_wide_search_is_never_worse_than_the_usual_one` counts strictly cheaper cases); that is how the missed edit above showed up.
- **Profiling the core from Python doesn't work with `--toggle-collect`** (2026-11-13): the pyo3 symbol names did not match and callgrind collected 83 instructions. Dump the core call's arguments from Python instead (a spy object replacing `nat._core` that writes `name;qubits;params` lines, plus the native gates and, for a device, `device <qubits> a-b ...`) and run them in a throwaway `core/examples/` driver (`compile_relabelling` or `device::compile_for_device`); move it out of `core/examples/` before `cargo clippy --all-targets`.
- **Byte-for-byte checks for network changes** (2026-11-13): the benchmark's three Hamiltonians are few; add held-out random Hamiltonians (other seeds, 6–24 qubits) to the outputs file, with CX and CZ.
- **Measuring a quality idea before building it** (2026-11-13): a pytket pass applied to QCC's output is a quick upper bound, but on a chip check the result's two-qubit gates are still on coupled pairs.
- **CI.** CI may run a newer Rust and clippy than the run's machine (a 1.99 clippy lint stopped a simulator release built with 1.97), so update to the latest stable before releasing. While master's version equals the one on PyPI, CI's `pip install --find-links dist` can pick the published wheel instead of the local build.
- **A keep rule that falls back can turn an improvement into a regression** (2026-11-15). "Keep A if it dominates B, else B" is not monotone: improving B in one count can make A stop dominating it, and the rule then picks a B with more CXs than the old A. When a new pass improves candidates that an existing rule compares, judge that rule on the costs without the new pass (`plain_cost` in `pipeline`).
- **Prototype a new rewrite on QCC's output in Python first** (2026-11-15): a numpy/Qiskit script that applies it to the final circuits and recompiles them (all-to-all and on the chip, checking couplings) measured the opportunity in 20 minutes (est. CNOTs saved per circuit), before any Rust; the Rust version then went through four placements in the pipeline, each measured over all 80 benchmark compilations with a throwaway env var (removed before timing).
- **Paired timing with two venvs** (2026-11-15): two long-lived worker processes (old wheel, new editable build) reading circuit keys on stdin, the driver alternating which goes first, medians of 9 (5 on the chip). Avoids process start-up and import noise.
- **Judging a pruning rule offline** (2026-11-17): trace every candidate's cost after each round (a throwaway thread-local tag naming the candidate, `eprintln!` under an env var) over many compilations once, then replay every rule and margin on the traces in Python. One 2-minute run answered what would have taken a build and a benchmark per rule.
- **Profiling one chip candidate without Python** (2026-11-17): dump the routed circuit of candidate 0 from `compile_for_device` (env var, throwaway) and run it through `qcc::passes::compile` (first pass plus `until_fixed_point`, as `on_the_same_pairs` does) in a throwaway `core/examples/` binary; it reproduced the candidate's result exactly (1231 / 4208).
- **Mutation-check new equivalence tests** (2026-11-17): reversing the order of a gate's Clifford steps went undetected because every gate's steps equal their reverse up to phase; dropping the inverses was caught. Pick a mutation that can matter.
