# qcc: the Rust core of qevedo-compiler

This crate does the compiler's heavy lifting. The Python package
(`../python`) parses OpenQASM, describes devices and calls into it through
PyO3 (`../python/src/lib.rs`); `../wasm` compiles it for the browser
(the decomposition playground on qevedo.com).

```bash
cargo test --release                              # all tests
cargo run --release --example multicore -- 100 10 # compare multi-core mappers
```

## Map of the code

| module | what it does |
|---|---|
| `ir.rs` | the instruction list every pass works on |
| `passes.rs` | the default pipeline: lowering, single-qubit merging, commutative cancellation, two-qubit block resynthesis |
| `lower.rs` | lowering one gate to the device's native gates (definitions or synthesis, whichever is cheaper) |
| `gates.rs` | matrices of the standard gates |
| `linalg.rs` | small fixed-size complex matrices (2×2, 4×4) |
| `one_qubit.rs` | optimal single-qubit synthesis (Euler angles) per native basis |
| `two_qubit.rs` | two-qubit synthesis: KAK decomposition, templates, moving rotations across gates |
| `numeric.rs` | numerically found decompositions for iSWAP-like native gates |
| `decompose.rs` | decompose one gate or unitary and check it (used by the playground) |
| `multicore/` | mapping circuits onto multi-core architectures: HQA and baselines (see `multicore/mod.rs` and `../docs/MULTICORE.md`) |
| `rng.rs` | a small seeded random generator |
| `testing.rs` | helpers for the unit tests |

The Python passes in `../python/qevedo/compiler/passes/` are a readable
reference implementation of the same pipeline; `../docs/ARCHITECTURE.md`
describes both.
