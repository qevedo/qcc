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

1. **DecomposeToNative** — lowers `h`, `u`, `u1`–`u3`, `p`, `x`, `y`, `z`, `s`, `t` into `{rz, sx, x}`; preserves `cx`, `measure`, `barrier`.
2. **CancelAdjacentInverses** — removes back-to-back self-inverses (`x`, `cx`, `sx`, opposite `rz`).

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

Reference profiles live in `python/profiles/`.

## Design choices

- **Python first** for IR and passes; C++ reserved for hot paths later (P2.7).
- **Frontend on `openqasm` 3.x**: parsing and semantic checks come from the `openqasm` package; `qcc/frontend/qasm.py` only lowers the checked tree to IR, and `qcc/emit/qasm.py` builds an `openqasm` tree and prints it, so output is always valid OpenQASM 2 or 3.
- **DeviceSpec drives routing**, not hard-coded IBM/Qiskit backend objects — this is the main differentiator vs UCC.
