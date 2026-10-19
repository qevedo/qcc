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

The default pipeline is: elide SWAPs, decompose, then two rounds of [phase
folding → merge single-qubit runs → commutative cancellation → two-qubit
block resynthesis], then moving gates through symmetric blocks and a final
merge. Each optimization can expose work for the others. For a device whose
two-qubit gate is CX-like but not CX (CZ, ECR, …), the Rust core also
compiles for CX first and then for the device, and keeps the cheaper result:
phase folding sees parities through CNOTs but not through the single-qubit
gates that make a CZ one (on the benchmark's adders, 40 and 92 CZs instead
of 49 and 113).

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
4. **MergeSingleQubitGates**: replaces each run of single-qubit gates on a
   qubit by its shortest equivalent in the native basis.
5. **CommutativeCancellation**: each gate looks back along its qubits,
   through gates it commutes with (Z-type gates on a CX control, X-type
   gates on a CX target, diagonal gates with each other), for an inverse to
   cancel or a rotation to merge with: `rz(a) q0; cx q0, q1; rz(b) q0`
   becomes `rz(a + b) q0; cx q0, q1`. (`CancelAdjacentInverses` is the
   simpler adjacent-only version.)
6. **ResynthesizeTwoQubitBlocks**: collects maximal runs of gates that only
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
7. **Symmetric blocks** (`core/src/symmetric_blocks.rs`, Rust only): the same
   symmetry lets a factor V ⊗ V move through a whole block, from the
   single-qubit runs before it to those after it. Blocks are synthesized
   without seeing their neighbours, so this pass picks V for each symmetric
   block to make the four runs around it cheapest. On Heisenberg chains
   (XX + YY + ZZ blocks) it saves 8–9% of the gates, and with the
   symmetric outer gates above it turned QCC's output for CZ from a fifth
   longer than Qiskit's to 7–13% shorter.
On the compilation benchmark (`benchmarks/compilation`: 21 circuits against
Qiskit's level 3, pytket and UCC), QCC uses no more two-qubit gates and no
more gates in all than Qiskit on any circuit, and fewer on most; phase
folding is most of the difference (19% fewer CNOTs on adders).

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

### Multi-core mapping

`core/src/multicore/` maps circuits onto modular architectures (several
cores linked by teleportation): Hungarian Qubit Assignment, FGP-rOEE and a
naive baseline, exposed as `qevedo.compiler.multicore`. It works on the
circuit's two-qubit gates grouped into timeslices and is separate from the
lowering pipeline. See [MULTICORE.md](MULTICORE.md).

### Planned passes (Phase 3)

- Initial layout (logical → physical qubit map)
- SABRE-style routing with SWAP insertion
- Technology-specific constraints from `DeviceSpec.gate_rules`

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
