# QCC (Python)

Technology-agnostic quantum circuit compiler. This is the active implementation; the C++ tree is reserved for future performance work.

## Install

```bash
cd qcc/python
pip install -e ".[dev]"
```

`openqasm` 3.x is installed from PyPI. To develop both together, also install
the sibling checkout: `pip install -e ../../openqasm`.

## Usage

```bash
qcc compile examples/bell.qasm -o /tmp/out.qasm --device profiles/all_to_all.yaml
qcc compile program.qasm -o /tmp/out.qasm --qasm-version 3 -I lib/
```

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
