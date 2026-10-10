# QCC: the Quantum Compiler Collection

A technology-agnostic quantum circuit compiler, part of the Qevedo packages:
it installs as `qevedo-compiler` and imports as `qevedo.compiler`.

## Install

```bash
pip install qevedo-compiler
```

For development, from a checkout (this builds the Rust core, so it needs a
Rust toolchain, from https://rustup.rs):

```bash
cd qcc/python
pip install -e ".[dev]"            # or: maturin develop --release
```

The Rust core lives in `../core` (`cargo test --release` there runs its own
tests); `src/` holds its Python bindings.

`openqasm` 3.x is installed from PyPI. To develop both together, also install
the sibling checkout: `pip install -e ../../openqasm`.

## Usage

```bash
qcc compile examples/bell.qasm -o /tmp/out.qasm --device qevedo/compiler/profiles/all_to_all.yaml
qcc compile program.qasm -o /tmp/out.qasm --qasm-version 3 -I lib/
```

From Python:

```python
from qevedo.compiler import compile_source, emit_qasm

circuit = compile_source(open("examples/bell.qasm").read())
print(emit_qasm(circuit))
```

The reference device profiles ship with the package, in
`qevedo/compiler/profiles/`.

Input can be OpenQASM 2 or 3. Programs are checked by `openqasm.analyze`
before lowering, so errors are reported with their location. The frontend
supports qubit and bit registers, gate calls with constant parameters
(broadcast over registers), measurements, resets, barriers, `let` aliases of
qubits and `const` values. Control flow, gate modifiers, classical variables,
timing and hardware qubits are rejected with an error until the IR supports
them. Gates defined by the program (or its own include files) are expanded
into their bodies; gates from `qelib1.inc` and `stdgates.inc` stay named so
the compiler can lower them optimally.

## Lowering to native gates

`compile_source` and `compile_file` rewrite every gate into the device
profile's `native_gates`, with the fewest gates:

- single-qubit gates through an Euler decomposition into the device's
  single-qubit basis (`rz`+`sx`, `rz`+`rx`, `rz`+`ry`, `rx`+`ry` or `u3`);
- two-qubit gates through a KAK decomposition into the device's two-qubit
  gate (`cx`, `cz`, `cy`, `ch`, `ecr`, `rxx`, `ryy`, `rzz` or `rzx`), using
  the minimal number of two-qubit gates (at most 3; on a `cx` device, 1 for
  `cz`, 2 for `crz` and 3 for `swap`);
- standard gates with a shorter textbook definition, and gates on three or
  more qubits, through their definitions.

After lowering, the circuit is optimized: single-qubit runs are merged,
gates cancel or merge across gates they commute with, and every run of gates
on one pair of qubits is resynthesized as a single two-qubit unitary when
that is cheaper. On a 12-qubit QFT and Qiskit's random circuits this gives the
same CX counts as Qiskit's `optimization_level=3` and equal or lower total
gate counts (see `docs/ARCHITECTURE.md`).

```python
from qevedo.compiler import CompileOptions, compile_source
from qevedo.compiler.device import DeviceSpec

ion_trap = DeviceSpec(name="ions", topology="all_to_all", qubits=4,
                      native_gates={"rx", "ry", "rxx"})
circuit = compile_source(open("program.qasm").read(), CompileOptions(device=ion_trap))
print(circuit.gate_counts())
```

## Tests

```bash
pytest
```
