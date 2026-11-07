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

1. **Quantum Volume speed, all-to-all.** QV-32 takes 243 ms against Qiskit's
   49 ms (2026-11-07). Serial passes after synthesis dominate: every
   instruction owns a `String` and three `Vec`s, so each copy or relabel
   allocates four times, and malloc/free are about a fifth of the
   instructions. The Python wrapper is 17–31 ms of it (profiles 2026-11-01
   and 2026-11-03).
2. **Hamiltonian compile time.** 0.40–0.45 s at 32 qubits against Qiskit's
   15 ms, after 0.4.7 added four trials. Ideas: stop a trial's pipeline once
   it can no longer beat the direct compilation; cache the templates' KAK
   per Weyl point.
3. **CZ rounds.** The fixed-point rounds (0.4.4) still cost most of the CZ
   time: the 16-qubit Hamiltonian on the chip takes 1.6 s against Qiskit's
   0.15 s. Each round's through-CX candidate leaves 100–600 more single-qubit
   gates for the next round to fold. The untried fix is a one-pass CX→CZ
   translation that's good enough to end the rounds sooner. See per-CZ H-side
   orientation under **Didn't work**: it's the closest lead.
4. **Extra two-qubit gates on the chip.** QV-8 uses 111 against Qiskit's 108.
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

**Speed (`lower.rs`, `parallel.rs`)**
- **A thread for every compilation** (2026-11-06). 8-qubit GHZ with CZ got 1.2× slower (1.5 → 1.7 ms): a thread costs a few tenths of a millisecond. There's now a 64-instruction threshold, and no thread without a second candidate.
- **Reading a profile as "no single hotspot"** (2026-10-30). The QV-32 profile looked spread out (resynthesis 44%). A callgrind look on 2026-11-01 found the serial passes were the bottleneck and the block syntheses on 4 cores were not. Profile the main thread, not totals, and note that under valgrind threads run one at a time.

## Considered, not tried

- **Python wrapper on Quantum Volume** (2026-11-01, 2026-11-03, 2026-11-06). About 17–31 ms of QV-32's time: `PassManager.run` copies the circuit and `NativePipeline` builds about 9,000 dataclasses (73,330 `__init__` calls over 5 runs). Deferred each time as a small share against the Rust core's 95–105 ms. It's still a candidate, behind **Next** item 1.
- **A new synthesis method instead of exact speed-ups** (2026-11-01). Speed-ups were chosen because they can be verified byte for byte.
- **LightSABRE (arXiv:2409.08368) and CZ-native single-qubit reduction** (2026-11-06). Searched; no decision recorded.
- **pytket `ThreeQubitSquash`** (2026-10-30). See **Next** item 6.

## Method lessons

- **Timing on the VM.**
  - Compare builds only side by side, interleaved, in the same run, with medians of 9–11.
  - Short runs (a few ms) need hundreds of calls.
  - With every core busy, any removed work counts, not just the critical path.
- **Tests at tolerance edges.** A count test exactly at the tolerance (θ − π/2 = 1e-10) can go either way under rounding in both methods. The 2026-11-01 run changed the test's nudges to either side of the edge. Never weaken a test to pass; if an edge is genuinely ambiguous, say so in the log.
- **Outputs.** For exact speed-ups, compare all 112 benchmark compilations byte for byte before and after.
- **CI.** CI may run a newer Rust and clippy than the run's machine (a 1.99 clippy lint stopped a simulator release built with 1.97), so update to the latest stable before releasing. While master's version equals the one on PyPI, CI's `pip install --find-links dist` can pick the published wheel instead of the local build.
