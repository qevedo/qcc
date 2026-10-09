# QCC: the Quantum Compiler Collection

A technology-agnostic quantum circuit compiler, part of the Qevedo packages:
it installs as `qevedo-compiler` and imports as `qevedo.compiler`.

## Install

```bash
pip install qevedo-compiler
```

For development, from a checkout:

```bash
cd qcc/python
pip install -e ".[dev]"
```

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
them. Gate definitions are not inlined: gates stay named instructions for the
decomposition pass.

## Tests

```bash
pytest
```
