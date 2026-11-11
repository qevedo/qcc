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

1. **Quantum Volume speed, all-to-all, continued.** QV-32 is 89 ms against
   Qiskit's 50 ms (2026-11-11, was 243). Two parts left:
   - The pyo3 crossing: the core call takes 70–75 ms from Python against
     59 ms in a Rust binary. pyo3 builds a list, a tuple and a string per
     instruction each way; interned name strings or flat arrays would cut
     it.
   - Synthesis: `TwoQubitSynth::dress`, `Shape::of`, `M4` products, `hypot`
     and `sincos` top the core's profile now. Allocation is down to 10% of
     the core's instructions (`phase_folding` groups, `dress`,
     `candidate_angles` remain).
2. **CZ device time.** The 16-qubit Hamiltonian on the chip takes 1.42 s
   (Qiskit 0.15 s). Ideas: a one-pass CX→CZ translation good enough to end
   the rounds sooner (see per-CZ H-side orientation under **Didn't work**);
   skip the second routing's rounds once its first round is clearly worse.
3. **Hamiltonian compile time.** 0.33 s at 32 qubits against Qiskit's
   17 ms. Ideas: stop a trial's pipeline once it can no longer beat the
   direct compilation; cache the templates' KAK per Weyl point.
4. **Routing quality beyond more trials.** QV-8 on the chip still uses 111
   CNOTs (Qiskit 108 at seed 1, but 108 or 111 depending on its seed);
   no search width tried here (up to 128 layouts or 32 routings each)
   reached 108. Untried from LightSABRE (Zou et al., arXiv:2409.08368):
   its release valve and relative scoring; routing whole two-qubit blocks
   as one node (tried as a lookahead change only, see **Didn't work**).
5. **The adder's gathered CNOTs** (37/85/181 at 8/16/32 qubits, against
   40/92/196 kept). They cost 108/252/540 gates against 105/245/525, so the
   keep rule declines them. Win them back by choosing the outer gates of the
   SWAP-elided blocks with their neighbours in view.
6. **Untested lead from 2026-10-30.** pytket's `ThreeQubitSquash` (with
   swaps) run on QCC's Clifford+T-32 output reached 306 CNOTs, against QCC's
   310 and pytket's own pipeline's 310. Three-qubit resynthesis was never
   tried in QCC.

## Worked

Newest first. Details are in `improvement-log.md` under the same date.

| Date | Version | Change | Effect |
|---|---|---|---|
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

**Gathering and single-qubit trades (`commuting_blocks.rs`, `symmetric_blocks.rs`)**
- **Moving a blocker out of the whole surrounding runs** (2026-10-30). The adder only went 40 → 39, because a √X earlier in a run stops the blocker. Cutting inside the runs and carrying third-qubit gates works.
- **Gathering without the keep rule** (2026-10-30). The adder reached 37/85/181 CNOTs but more gates in all. The pipeline now runs with and without gathering and keeps the cheaper.
- **CX trades before routing, or in the CX stage of CZ compilation** (2026-10-30). Chip CZ Clifford+T-8 cost 100 more gates, because phase folding lost its Z rotations. Trades run only in the last pipeline (`last`).
- **A guard keeping diagonal runs diagonal** (2026-10-30). No gain (CZ Clifford+T-8 still 236/965). Running the trades last was the fix.
- **Accepted weaknesses after 0.4.2** (2026-10-30). Some chip rows got slightly worse in one count from routing variance: CX Clifford+T-16, CZ adder-16, CZ Hamiltonian-8, CZ QFT-16. Remedy, not yet tried: route the logical circuit with and without gathering and keep the cheaper.

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

**Speed (`lower.rs`, `parallel.rs`, `ir.rs`, Python wrapper)**
- **Ruled out while profiling the wrapper** (2026-11-11). Building the raw tuples (5.5 ms) and the output `Instruction` dataclasses (5 ms) on QV-32 are not the wrapper's cost; the garbage collector was (15–20 ms, see **Worked** 0.4.9). The 2026-11-01/03/06 entries that blamed the dataclasses were misread.
- **A thread for every compilation** (2026-11-06). 8-qubit GHZ with CZ got 1.2× slower (1.5 → 1.7 ms): a thread costs a few tenths of a millisecond. There's now a 64-instruction threshold, and no thread without a second candidate.
- **Reading a profile as "no single hotspot"** (2026-10-30). The QV-32 profile looked spread out (resynthesis 44%). A callgrind look on 2026-11-01 found the serial passes were the bottleneck and the block syntheses on 4 cores were not. Profile the main thread, not totals, and note that under valgrind threads run one at a time.

## Considered, not tried

- **Python wrapper on Quantum Volume** (2026-11-01, 2026-11-03, 2026-11-06). Done on 2026-11-11 in another way: most of it was the garbage collector, now paused (`paused_gc`). Left: `PassManager.run`'s circuit copy (3.8 ms on QV-32), building the tuples and dataclasses (about 5 ms each way), and pyo3's own conversion (**Next** item 1).
- **A global allocator (mimalloc, jemalloc) in the Python extension** (2026-11-11). Would cut the remaining 10% allocation cost without code changes, but adds a C dependency to every wheel target of the Release workflow and cannot serve the WebAssembly build; removing the allocations was exact and portable. Retry if allocation is again a large share after the remaining sources (phase folding's groups, `dress`, `candidate_angles`) are gone.
- **`smallvec` / `rustc-hash` crates** (2026-11-11). The core has a single dependency (`num-complex`); the two small types it needed (`Inline`, `FxHasher`) are about 250 lines with tests in `ir.rs` and `hashing.rs`.
- **A new synthesis method instead of exact speed-ups** (2026-11-01). Speed-ups were chosen because they can be verified byte for byte.
- **LightSABRE (Zou, Treinish, Hartman, Ivrii and Lishman, arXiv:2409.08368)** (2026-11-06, 2026-11-09). Read for 0.4.8: QCC already has its release valve and uses SABRE's nondeterminism (trials); 0.4.8 took its lesson that more trials pay once routing is cheap. Its relative scoring would change float ties, so the outputs could not be checked byte for byte; left for a run that aims at routing quality (Next item 4).
- **CZ-native single-qubit reduction** (2026-11-06). Searched; no decision recorded.
- **pytket `ThreeQubitSquash`** (2026-10-30). See **Next** item 6.

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
- **CI.** CI may run a newer Rust and clippy than the run's machine (a 1.99 clippy lint stopped a simulator release built with 1.97), so update to the latest stable before releasing. While master's version equals the one on PyPI, CI's `pip install --find-links dist` can pick the published wheel instead of the local build.
