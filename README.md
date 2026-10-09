# qcc

[![Python tests][qcc-ti]][qcc-tu]

Quantum Compiler Collection — a technology-agnostic quantum circuit compiler.
The Python package is published as `qevedo-compiler` and imports as
`qevedo.compiler`; the command is `qcc`.

## Status

The **Python compiler** in `python/` is the active implementation (v0.1):

- QCC IR (`Circuit`, `Instruction`, `DeviceSpec`)
- OpenQASM 2 and 3 frontend and emitter, built on [`openqasm`](https://pypi.org/project/openqasm) 3.x (sibling repo [`openqasm`](../openqasm))
- Decomposition to `{rz, sx, x, cx}` and basic optimization
- CLI: `qcc compile input.qasm -o out.qasm --device qevedo/compiler/profiles/grid_2x4.yaml [--qasm-version 3]`

The C++ tree is a placeholder for future performance work.

## Quick start

```bash
cd python
pip install -e ".[dev]"            # add -e ../../openqasm to use the sibling checkout
pytest
qcc compile examples/bell.qasm -o /tmp/bell_out.qasm --stats
```

See [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) and [python/README.md](python/README.md).

## Releases

Pushing a tag `vX.Y.Z` that matches `python/pyproject.toml` runs the tests,
builds the package and publishes it to PyPI through trusted publishing
(`.github/workflows/release.yml`).

## License

MIT. Copyright [Qevedo](https://qevedo.com)

[qcc-ti]: https://github.com/qevedo/qcc/actions/workflows/python.yml/badge.svg?branch=master
[qcc-tu]: https://github.com/qevedo/qcc/actions/workflows/python.yml
