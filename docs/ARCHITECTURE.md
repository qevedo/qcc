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

The default pipeline is decompose → merge single-qubit runs → cancel inverses →
merge single-qubit runs again.

1. **DecomposeToNative**: rewrites every gate that is not in the device's
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
2. **MergeSingleQubitGates**: replaces each run of single-qubit gates on a
   qubit by its shortest equivalent in the native basis.
3. **CancelAdjacentInverses**: removes pairs of gates that undo each other
   with no other gate on their qubits in between (`x x`, `cx cx`, `s sdg`,
   `rz(a) rz(-a)`, …; `sx sx` is an `x` and is left for the merge pass).

The synthesis lives in `qevedo/compiler/synthesis/` (`gates.py`: matrices of
the standard gates; `one_qubit.py`: Euler decompositions; `two_qubit.py`: KAK
and templates). Not supported yet as native two-qubit gates: `iswap`,
`sqrt_iswap` and other gates that are not CX-like or Ising rotations.

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

- **Python first** for IR and passes; C++ reserved for hot paths later (P2.7).
- **Frontend on `openqasm` 3.x**: parsing and semantic checks come from the `openqasm` package; `qevedo/compiler/frontend/qasm.py` only lowers the checked tree to IR, and `qevedo/compiler/emit/qasm.py` builds an `openqasm` tree and prints it, so output is always valid OpenQASM 2 or 3.
- **DeviceSpec drives routing**, not hard-coded IBM/Qiskit backend objects — this is the main differentiator vs UCC.
