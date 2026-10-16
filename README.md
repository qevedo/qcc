# QCC

[![Python tests][qcc-ti]][qcc-tu]

Quantum Compiler Collection — a technology-agnostic quantum circuit compiler.
The Python package is published as `qevedo-compiler` and imports as
`qevedo.compiler`; the command is `qcc`.

Try it in the browser: the [gate decomposition](https://qevedo.com/decompose/)
and [multi-core mapping](https://qevedo.com/multicore/) playgrounds run the
Rust core compiled to WebAssembly.

## Status

`qevedo-compiler` 0.2 on PyPI:

- QCC IR (`Circuit`, `Instruction`, `DeviceSpec`)
- OpenQASM 2 and 3 frontend and emitter, built on [`openqasm`](https://pypi.org/project/openqasm) 3.x ([repository](https://github.com/qevedo/openqasm))
- Lowering to any device's native gates with the fewest gates (Euler and KAK synthesis, optimal two-qubit gate counts) and peephole optimization
- CLI: `qcc compile input.qasm -o out.qasm --device qevedo/compiler/profiles/grid_2x4.yaml [--qasm-version 3]`
- Multi-core (modular) quantum computers ([docs/MULTICORE.md](docs/MULTICORE.md)): mapping with Hungarian Qubit Assignment (HQA), a latency-aware HQA, FGP-rOEE and a naive mapper on any connected links between cores; routing inside each core with SABRE, optionally through fixed communication qubits; a timing model of links with a limited supply of EPR pairs; and benchmarks reproducing the published comparison and comparing with TeleSABRE
- A preprint draft on the multi-core work, generated from the benchmark results ([paper/](paper/))

The lowering and optimization passes run in a **Rust core** (`core/`, the
`qcc` crate), loaded by the Python package as `qevedo.compiler._core` through
PyO3 bindings (`python/src/`). The Python passes stay as the reference
implementation, and the tests check one against the other.

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
