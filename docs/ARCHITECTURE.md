# QCC IR

QCC uses a small, owned intermediate representation designed to be independent of Qiskit, Cirq, or PyTKET.

## Core types

| Type | Role |
|------|------|
| `Circuit` | Container: register sizes + ordered `Instruction` list |
| `Instruction` | Gate/operation name, qubit operands, optional parameters and classical bits |
| `Qubit` / `Clbit` | Named register + index (OpenQASM-style) |
| `DeviceSpec` | Target topology, native gateset, scheduling hints |

## Pipeline (v0.1)

```
OpenQASM 2/3  →  openqasm.parse + analyze  →  frontend  →  IR  →  passes  →  emitter  →  OpenQASM 2/3
```

### Implemented passes

The default pipeline is: elide SWAPs, decompose, [rotation merging → phase
folding → merge single-qubit runs → commutative cancellation → gathering
blocks through commuting gates → two-qubit block resynthesis], the same again without resynthesis, then moving gates through
symmetric blocks and a final merge. Each optimization can expose work for
the others; a second resynthesis never removed a gate on the benchmarks, and
resynthesis is most of the time. The Rust core synthesizes a pass's blocks
on all cores when there are many (`Lowering::prefetch_two_qubit`), and runs
work that reads only the same input side by side (`core/src/parallel.rs`):
the direct and through-CX compilations below, the Pauli networks beside
them, SABRE's layout trials, and the routing of each logical candidate on a
chip. The results come back in the order a loop would produce them, and the
first of the cheapest is kept, so the output does not depend on the
threads. All the pipelines of one compilation share their two-qubit
syntheses (`SynthesisCache` in `core/src/lower.rs`): the rounds below see
mostly the same blocks again, and on the 16-qubit Hamiltonian on the
heavy-hex chip with CZ two thirds of the syntheses had been done before.
A synthesis depends only on the unitary, its objective and the native
gates, so sharing them changes no output. For a device whose
two-qubit gate is CX-like but not CX (CZ, ECR, …), the Rust core also
compiles for CX first and then for the device, and keeps the cheaper result:
phase folding sees parities through CNOTs but not through the single-qubit
gates that make a CZ one (on the benchmark's adders, 40 and 92 CZs instead
of 49 and 113). For such devices the pipeline is then run again on its own
result, without relabelling or new pairs, while that gets cheaper (at most
eight rounds; `until_fixed_point` in `core/src/passes.rs`): each translation
back from CX leaves single-qubit gates where the next reading through CX
finds more to fold. On the heavy-hex chip with CZ it took the 16-qubit
Hamiltonian from 1379 CZs and 5372 gates to 1292 and 4572.

1. **ElideSwaps** (`passes/permutation.py`; `core/src/permutation.rs`):
   removes SWAP gates by relabelling the qubits of everything after them, as
   Qiskit's `ElidePermutations` does. That leaves qubits' states on other
   wires at the end, which cannot be seen when each moved state is measured
   last; otherwise the pass only runs with `CompileOptions(relabel_outputs=True)`,
   and the result's `output_permutation` says where each state ends. It
   leaves circuits for devices that are not all-to-all alone, since
   relabelled gates may act on uncoupled pairs.
2. **DecomposeToNative**: rewrites every gate that is not in the device's
   `native_gates` with the fewest native gates, keeping the cheapest of the
   routes that apply (fewest two-qubit gates first, then fewest gates):
   - *Single-qubit gates*: a ZYZ Euler decomposition emitted in the device's
     single-qubit basis (`rz`+`sx`[+`x`], `rz`+`rx`, `rz`+`ry`, `rx`+`ry` or
     `u3`), which is optimal for each basis.
   - *Two-qubit gates*: a KAK decomposition. The Weyl-chamber point of the
     gate fixes the minimal number of native two-qubit gates (0–3 for CX-like
     gates `cx`, `cz`, `cy`, `ch`, `ecr`; one per non-zero coordinate for
     `rxx`, `ryy`, `rzz`, `rzx`). A template circuit with that point is
     dressed with single-qubit gates computed from the two KAK
     decompositions; template variants and rotations moved across the
     two-qubit gates (Z on a CX control, X on its target, Pauli products
     through Clifford gates) keep the single-qubit gate count low.
   - *Standard-library definitions*: `qelib1.inc`/`stdgates.inc` gates with a
     known short definition (`crz` is `rz, cx, rz, cx`), lowered recursively.
     Gates on three or more qubits (`ccx`, `cswap`, `rccx`, and the `c3x`
     family through their `qelib1.inc` bodies) take this route.
3. **PhaseFolding** (`passes/phase_folding.py`; `core/src/phase_folding.rs`):
   merges rotations that act on the same parity of the qubits, however far
   apart (Amy, Maslov and Mosca 2014; Nam et al. 2018). Each wire carries an
   affine parity of variables: CNOTs XOR parities, X flips their constant,
   and any other non-diagonal gate gives its qubits fresh variables. Diagonal
   rotations (Rz, phase, T, S, Z, Rzz) add phases on their parity, and those
   on the same parity merge into the last adjustable one.
   In the Rust core, **rotation merging** (`merge_rotations` in
   `core/src/pauli_network.rs`) goes first and further: each Z rotation,
   moved through every Clifford gate before it (Hadamards included, with a
   stabilizer tableau), is a rotation about a Pauli string, and rotations
   about the same string merge when every rotation between them commutes
   with it (Zhang and Chen, 2019).
4. **MergeSingleQubitGates**: replaces each run of single-qubit gates on a
   qubit by its shortest equivalent in the native basis.
5. **CommutativeCancellation**: each gate looks back along its qubits,
   through gates it commutes with (Z-type gates on a CX control, X-type
   gates on a CX target, diagonal gates with each other), for an inverse to
   cancel or a rotation to merge with: `rz(a) q0; cx q0, q1; rz(b) q0`
   becomes `rz(a + b) q0; cx q0, q1`. (`CancelAdjacentInverses` is the
   simpler adjacent-only version.)
6. **Gathering blocks** (`core/src/commuting_blocks.rs`, Rust only): two
   gates on the same pair that a gate on one of their qubits and a third
   separates fall in different blocks, even when that gate commutes with
   one of them: `cx a,b; cx a,c; cx b,a` is two blocks, but `cx a,c`
   commutes with `cx a,b`, and the two CXs on `{a, b}` are a SWAP and one CX.
   This pass moves such gates out of the way, as tket's Clifford
   simplification does for CXs (Sivarajah et al., arXiv:2003.10611; Fagan
   and Duncan, arXiv:1901.10114), but for any two-qubit gate and with an
   exact check: a gate moves only if it commutes with the product of the
   pair's gates it crosses (an 8×8 check, so single-qubit gates count
   together), taking along the gates on its third qubit that depend on it,
   and only where the gathered block needs fewer CX-like gates. A gathered
   block can trade single-qubit gates for CXs (the SWAP elided from it
   leaves generic single-qubit gates behind), so the pipeline finishes the
   circuit with and without gathering whenever gathering moved a gate, and
   keeps the gathered one only when it needs no more two-qubit gates and
   no more gates in all.
   **Carrying gates between blocks** (`carry_through` in the same file):
   where the gates between two blocks on a pair do not commute with them,
   part of a block may, once moved through the single-qubit gates between.
   In `cx a,b; h b; cx b,c; cx a,b`, `h b · cx a,b` is `cz a,b · h b`, and
   the CZ commutes with the control of `cx b,c`: it joins the last CX, and
   `cx·cz` on `{a, b}` is one CX, so the circuit needs two CXs instead of
   three. The pass carries the end of a block (from one of its last three
   two-qubit gates) forward into the next block on its pair, or the start
   of a block back into the previous one: a carried part `D` crosses a
   single-qubit gate `g` as `g D g†` and a blocker only if it commutes
   with `D` as it stands there (an 8×8 check), and is put down as its own
   gates between the inverse of the single-qubit gates it crossed and
   those gates again, for the resynthesis to merge. A move is made only
   when the two blocks need fewer CX-like gates. It runs only on a result
   that will not be optimized again (after routing, on the chip's qubits),
   after block resynthesis; the blocks it changes are synthesized again,
   and the result is kept only when it needs no more two-qubit gates and
   no more gates in all than without it. Gathering's choice above is made
   on the costs without carrying, so carrying cannot change it. For CZ-like
   devices the repeated rounds (`until_fixed_point`, above) first run
   without carrying until they stop shrinking, then with it, so that the result is never costlier than
   before. On the heavy-hex chip it took the 16-qubit Hamiltonian from
   1336 to 1321 CNOTs and, with CZ, from 1256 / 4454 to 1231 / 4208 CZs /
   gates.
7. **ResynthesizeTwoQubitBlocks**: collects maximal runs of gates that only
   touch one pair of qubits, multiplies each into its 4×4 unitary and
   replaces it by a fresh KAK synthesis when that has fewer two-qubit
   gates, or as many and fewer gates. Two SWAPs on the same pair vanish;
   four CXs that amount to two become two. A block's outer single-qubit
   gates merge with its neighbours' afterwards, so the synthesis here
   minimises the gates between its first and last two-qubit gates first
   (`Objective::InCircuit`): saving an outer gate in isolation usually saves
   nothing once merged, and on Quantum Volume that alone was 1% of the gates.
   Where qubits may be relabelled (as for ElideSwaps), the Rust core's
   resynthesis is also permutation-aware: a block U whose SWAP·U needs fewer
   CX (a CX each way is a SWAP and one CX) becomes SWAP·U, its SWAP elided,
   and everything after it is relabelled (`compile_relabelling` in
   `core/src/passes.rs`; the Python passes do not do this).
   At a symmetric Weyl point the KAK decomposition's outer gates are not
   unique: canonical(a, b, c) commutes with Rz ⊗ Rz when a = b, with
   Rx ⊗ Rx(±) when b = ±c, and with every V ⊗ V (or V ⊗ ZVZ) when
   a = b = ±c, as for XX + YY + ZZ. The Rust core tries the outer gates those
   symmetries allow and keeps the cheapest (`symmetric_variants` in
   `core/src/two_qubit.rs`).
8. **Symmetric blocks** (`core/src/symmetric_blocks.rs`, Rust only): the same
   symmetry lets a factor V ⊗ V move through a whole block, from the
   single-qubit runs before it to those after it. Blocks are synthesized
   without seeing their neighbours, so this pass picks V for each symmetric
   block to make the four runs around it cheapest. On Heisenberg chains
   (XX + YY + ZZ blocks) it saves 8–9% of the gates, and with the
   symmetric outer gates above it turned QCC's output for CZ from a fifth
   longer than Qiskit's to 7–13% shorter. A block at the point (a, 0, 0),
   a lone CX for one, commutes with an X rotation on either qubit alone in
   its canonical frame (a Z rotation on a CX's control, an X rotation on
   its target), so each of its qubits trades its own rotation, chosen to
   make the runs before and after it cheapest; nothing else chooses those
   angles. It changes runs of Z rotations into generic ones, which phase
   folding cannot merge, so it runs only on a circuit that will not be
   optimized again (not before routing, nor in the first half of compiling
   through CX for a CZ device). With gathering, it took the benchmark's
   random Clifford+T circuits from 293, 525 and 1074 gates to 278, 505 and
   1034.
9. **Pauli pairs** (`core/src/pauli_pairs.rs`, Rust only): where a Pauli P
   on one qubit commutes with a two-qubit gate G (X with a CX's target, Z
   with its control or with either qubit of a CZ; checked on G's matrix),
   G = P G P, so a P can be put just before G and another just after it,
   each merging into the run on its side. Done around several gates along
   a wire, this carries a Pauli from one run to another: an X that a
   block's synthesis left alone between two CXs it is the target of can
   travel to a run further on that absorbs it for free (in {rz, sx, x},
   X·Rz(b)·SX·Rz(a) = Rz(π − b)·SX·Rz(π + a)), or meet another X and
   cancel. Pairs on different wires are independent, so each wire is
   solved exactly: a dynamic programme over its chain of runs, with the
   Pauli at each gate as the state, picks the cheapest choices, and a wire
   changes only if that is strictly cheaper. Qiskit's
   `Optimize1qGatesSimpleCommutation` moves gates one commutation at a
   time into runs that exist; this pass crosses any number of gates and can
   create a run. It runs on final results, at the end of each pipeline
   (so gathering above is judged with it), and for CZ-like devices once
   after the repeated rounds, as inside them it changed which fixed point
   they reached. With it, gathering pays on the adder: 37 / 104, 85 / 244
   and 181 / 524 at 8, 16 and 32 qubits, against 40 / 105, 92 / 245 and
   196 / 525.
   Since 0.4.14 the last pass of a compilation also offers rotations: a
   whole axis commutes with a CX's target (X) and control (Z) and with a
   CZ (Z), so G = D† G D for any rotation D about it, and the chain's
   state at each gate can be such a rotation. The angles offered at a gate
   are those that would remove an outer rotation (in that axis) of one of
   the two nearest runs on either side that it reaches through gates of
   the same axis, with their shifts by π and −π/2 (the other forms of the
   {rz, sx} decomposition), besides the identity and the Pauli, so the
   result is never costlier than with Paulis alone. A Z rotation can thus
   leave the end of one generic run and join the start of another through
   several CX controls or CZs, and an X rotation along CX targets. It runs
   once, on the result no pipeline will see again (`with_rotation_pairs`
   in `passes.rs`): inside the pipeline it changed the circuits that
   gathering, carrying and the choice between candidates compare, and
   their keep rules then fell back to circuits with more two-qubit gates
   (the adder lost its gathered CNOTs). Run last, it leaves every choice
   as it was: no two-qubit gate changes, and no circuit gets more gates.
   Between two Z-axis gates a run's {rz, sx} count is a constant plus one
   term for each side (`Basis1q::split_count`), so the cheapest choice
   before is found once for all choices after it. On the benchmark's
   32-qubit outputs the pass takes 0.2–0.8 ms, and 3–6 ms on Quantum
   Volume's 9 000–11 000 gates (about 1.2–2.7× the Pauli pass).
On the compilation benchmark (`benchmarks/compilation`: 24 circuits against
Qiskit's level 3, pytket and UCC), QCC uses no more two-qubit gates than
Qiskit on any circuit, and fewer on most; phase folding is most of the
difference on arithmetic (with gathering and Pauli pairs, 25% fewer
CNOTs on adders), the Pauli network on
Hamiltonian simulation (up to 43%).

The synthesis lives in `qevedo/compiler/synthesis/` (`gates.py`: matrices of
the standard gates; `one_qubit.py`: Euler decompositions; `two_qubit.py`: KAK
and templates), mirrored in Rust in `core/src/`. The Rust core also accepts
`cp` as a native gate (an Ising rotation up to Z rotations) and any other
entangling gate (`iswap`, `sqrt_iswap`, `sqrt_swap`, …): `core/src/numeric.rs`
finds the fewest uses numerically (Levenberg–Marquardt on the Makhlin
invariants, then a polish against the target matrix that makes the result
exact). The Python reference has templates for CX-like and Ising gates only.

`core/src/decompose.rs` decomposes one gate or unitary and checks the result;
`wasm/` compiles it to WebAssembly for the decomposition playground at
https://qevedo.com/decompose/.

### Pauli networks

For circuits whose two-qubit gates are mostly CX ladders around few
rotations, as Pauli exponentials in Hamiltonian simulation and chemistry
are, the Rust core also synthesizes the circuit again from scratch
(`core/src/pauli_network.rs`) and keeps the result if it is cheaper:

1. The circuit, lowered to CX, is read as Pauli rotations: each Z rotation
   moved through the Clifford gates before it (a stabilizer tableau) is a
   rotation about some Pauli string, and the Clifford gates that remain
   form one Clifford at the end.
2. The rotations are built again greedily, as in Rustiq (Goubault de
   Brugière and Martiel, 2024) and pytket's GreedyPauliSimp: of the
   controlled-Pauli gates on two qubits of a rotation that can come next,
   the one that most lowers the weights of the rotations ahead (discounted
   by commuting layer) and of the Clifford left; a rotation is emitted as an
   Rz once it acts on one qubit. The local Cliffords that set up a CX are
   not undone after it (the rotations ahead are read in the new frame, with
   the same weights), each CX is placed whichever way round needs fewer of
   them, and basis changes use √X rather than H: together they take about
   half the network's single-qubit gates off.
3. The Clifford left is synthesized one qubit at a time, as in the greedy
   synthesis of Bravyi et al. (2021); where qubits may be relabelled, it
   may end as a permutation.
4. The greedy choice is blind between gates that do about as well, so eight
   syntheses run side by side (one thread each), seven of them breaking
   those ties with seeded noise, as pytket's GreedyPauliSimp does with
   `trials`. Four more charge ten times as much for each local Clifford a
   candidate CX needs (a tenth of a weight instead of a hundredth), which
   gives up a few CXs for fewer single-qubit gates. Each trial, after the
   pipeline, is judged against the direct compilation on its own, and the
   first of the cheapest that pays is kept, the same on any machine
   (`network_trials` and `paying_network` in `core/src/passes.rs`).

Scoring is most of the synthesis's work: nine candidate gates (σ on one
qubit controlling τ on the other) on up to 64 pairs for every gate placed.
What a candidate does to the letters of its two qubits is the same on any
pair, so each candidate is a table of what becomes of the 16 pairs of
letters, worked out once (`candidate_tables`), and the nine scores of a
pair are summed in one pass over the rotations ahead, each in the same
order as nine passes, so they are the same numbers (since 0.4.10; the
32-qubit Hamiltonian compiles in a third of the time it took before).

It runs only on circuits with three or more CXs per non-Clifford rotation
(at most 400 rotations) and on all-to-all devices, and its result is kept
only when it needs no more two-qubit gates and no more gates in all than
the direct compilation. On one Trotter step of random Hamiltonians QCC then
uses 47.1% fewer CXs than Qiskit's level 3 at 8 qubits, 20.4% at 16 and
7.2% at 32, and 39.7%, 18.7% and 15.5% fewer gates in all (101 / 286, 313 / 790 and
687 / 1511 CXs / gates).

### Multi-core mapping

`core/src/multicore/` maps circuits onto modular architectures (several
cores linked by teleportation): Hungarian Qubit Assignment, FGP-rOEE and a
naive baseline, exposed as `qevedo.compiler.multicore`. It works on the
circuit's two-qubit gates grouped into timeslices and is separate from the
lowering pipeline. See [MULTICORE.md](MULTICORE.md).

### Compiling for a device

For a device that is not all-to-all and lists its couplings
(`DeviceSpec.edges`), the circuit is placed and routed onto its physical
qubits (`core/src/device.rs`):

1. The default pipeline runs on the logical qubits first, where relabelling
   is free: fewer two-qubit gates leave the router less to do.
2. The layout: if the pairs the circuit's gates use fit the coupling graph,
   a layout that needs no SWAP (`core/src/routing/perfect.rs`, a subgraph
   search like Qiskit's VF2Layout); otherwise SABRE's forward–backward
   search over 32 random starts.
3. Routing with SABRE, where a SWAP right after a block on the same pair
   counts for what it costs once merged into the block, as in NASSC (Liu,
   Li and Zhou, HPCA 2022): after a CX one CNOT, after a generic block none.
   Each layout is routed 16 times with different tie-breaking, and the
   routing with the cheapest SWAPs kept; so is the one the first 4 of each
   layout would have given (SABRE's default, and QCC's before 0.4.8).
4. The default pipeline again, on the physical qubits, keeping every
   two-qubit gate on its pair; on a device whose couplings have a direction,
   gates the other way round are turned around with single-qubit gates
   found from their KAK decompositions. When step 3 kept two routings, both
   are optimized and the cheaper kept: the SWAPs' CNOTs are an estimate of
   what the optimized circuit needs, so the wider search alone could make
   a circuit worse.
5. For a circuit of Pauli rotations (where the all-to-all pipeline tries
   Pauli networks), one more candidate: from each initial layout of
   step 3, the rotations placed on the physical qubits and synthesized
   again with every CX that sets one up on a coupled pair, scored by each
   rotation's spanning tree on the chip (`synthesize_on` in
   `core/src/pauli_network.rs`), the Clifford left at the end routed
   after; the trial with fewest CXs and SWAPs after a light cleanup goes
   through step 4 beside the others (`networks_on_chip` in
   `core/src/device.rs`). The 16-qubit Hamiltonian on the heavy-hex chip:
   1321 → 1078 CNOTs (0.4.15).

The result is a circuit on a register `q` of the device's qubits;
`Circuit.initial_layout` and `Circuit.output_permutation` say on which
physical qubit each logical qubit starts and ends. A circuit whose two-qubit
gates already act on coupled pairs keeps its layout. On a 19-qubit
heavy-hex chip (`benchmarks/compilation/device.py`) this uses 12–14% fewer
two-qubit gates than Qiskit's level 3 (0.880× with CX, 0.862× with CZ), and
19–21% fewer gates in all.

Still to come: technology-specific constraints from `DeviceSpec.gate_rules`,
and error rates in the layout.

## DeviceSpec schema

```yaml
name: grid_2x4
topology: grid          # coupling_graph | grid | linear_chain | all_to_all | ion_trap_zones
qubits: 8
edges: [[0, 1], [1, 2]] # undirected unless directed_edges: true
native_gates: [rz, sx, x, cx]
directed_edges: false
max_parallel_2q: 2
gate_rules:
  cx: { on: edge }
```

Reference profiles live in `python/qevedo/compiler/profiles/` and ship with the package.

## Design choices

- **Rust core, Python reference.** The default pipeline (`PassManager.default`)
  is one pass, `NativePipeline`, that hands the instruction list to the Rust
  core (`core/`, crate `qcc`, imported as `qevedo.compiler._core`). The core
  ports the Python passes one to one: the same gate matrices, Euler and KAK
  synthesis, templates, rotation moves, definitions and pipeline. The Python
  passes (`PassManager.reference`) remain the readable reference;
  `tests/test_native.py` checks that both give the same two-qubit counts gate
  by gate and similar totals on whole circuits. The two can differ by a gate
  or two per gate where the rotation search lands in different local optima.
  The Rust pipeline is 40–80× faster (31 ms against 1.3 s on a 10-qubit,
  855-gate random circuit), on par with Qiskit's transpiler.
- **Frontend on `openqasm` 3.x**: parsing and semantic checks come from the `openqasm` package; `qevedo/compiler/frontend/qasm.py` only lowers the checked tree to IR, and `qevedo/compiler/emit/qasm.py` builds an `openqasm` tree and prints it, so output is always valid OpenQASM 2 or 3.
- **DeviceSpec drives routing**, not hard-coded IBM/Qiskit backend objects — this is the main differentiator vs UCC.
