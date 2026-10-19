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

After lowering, the circuit is optimized: rotations acting on the same
parity of the qubits are merged however far apart they are (phase folding),
single-qubit runs are merged, gates cancel or merge across gates they
commute with, and every run of gates on one pair of qubits is resynthesized
as a single two-qubit unitary when that is cheaper, choosing among equally
short syntheses the one with the fewest gates inside the block. Two-qubit
interactions whose removal costs less than 10⁻⁹ in average gate infidelity
are dropped, as Qiskit's default does. On a benchmark of 21 circuits (QFT,
QAOA, Quantum Volume, GHZ, Heisenberg, adders, Clifford+T) this uses 0.968×
the two-qubit gates of Qiskit's `optimization_level=3` and 0.901× its gates
in all, and no more of either on any circuit
([benchmarks/compilation](https://github.com/qevedo/qcc/tree/master/benchmarks/compilation)).

SWAP gates are removed by relabelling the qubits after them, and blocks that
are a SWAP away from fewer CX are synthesized that way, whenever the moved
qubits are measured last, so that the output is equivalent. With
`CompileOptions(relabel_outputs=True)` (`qcc compile --relabel-outputs`)
QCC does so for any circuit, as Qiskit does, and the result's
`output_permutation` says on which wire each qubit's state ends. On a device
that is not all-to-all, see below.

```python
from qevedo.compiler import CompileOptions, compile_source
from qevedo.compiler.device import DeviceSpec

ion_trap = DeviceSpec(name="ions", topology="all_to_all", qubits=4,
                      native_gates={"rx", "ry", "rxx"})
circuit = compile_source(open("program.qasm").read(), CompileOptions(device=ion_trap))
print(circuit.gate_counts())
```

## Compiling for a device

For a device whose qubits are not all coupled, give its couplings: QCC
optimizes the circuit, places it (with no SWAP at all when the circuit's
pairs fit the chip), routes it with SABRE, counting each SWAP for what it
costs once merged with the gates next to it, and optimizes it again on the
physical qubits, turning gates around on devices whose couplings have a
direction.

```python
from qevedo.compiler import CompileOptions, compile_source
from qevedo.compiler.device import DeviceSpec

chip = DeviceSpec(name="line", topology="linear_chain", qubits=5,
                  edges=[(0, 1), (1, 2), (2, 3), (3, 4)],
                  native_gates={"rz", "sx", "x", "cz"})
circuit = compile_source(open("program.qasm").read(), CompileOptions(device=chip))
print(circuit.initial_layout, circuit.output_permutation)  # where each qubit starts and ends
```

The result is on the device's physical qubits (a register `q`), and
`qcc compile --device profile.yaml` prints both layouts. On a 19-qubit
heavy-hex chip, QCC alone uses 5–6% fewer two-qubit gates than Qiskit's
`optimization_level=3` and 12–17% fewer gates in all
([benchmarks/compilation](https://github.com/qevedo/qcc/tree/master/benchmarks/compilation)).

## Inside Qiskit

From version 0.2.1, with Qiskit installed (`pip install qevedo-compiler[qiskit]`), QCC is also
an optimization stage for Qiskit's transpiler: Qiskit lays out, routes and
translates the circuit, and QCC optimizes it.

```python
from qiskit import transpile

compiled = transpile(circuit, backend, optimization_method="qcc")
```

On a 19-qubit heavy-hex chip with CX gates, this uses 3% fewer two-qubit
gates than Qiskit's level 3 and 14% fewer gates in all (61% on a QFT); with
CZ gates, 9% fewer gates in all; see
[benchmarks/compilation](https://github.com/qevedo/qcc/tree/master/benchmarks/compilation).

## Multi-core architectures

`qevedo.compiler.multicore` compiles circuits for modular quantum computers
made of several cores joined by quantum links, in Rust:

- **mapping**: before each timeslice, which core every qubit is in, with as
  few teleportations as possible: Hungarian Qubit Assignment (HQA, Escofet
  et al.), FGP-rOEE and a naive mapper, on cores linked all to all, in a
  line, a ring, a grid or any connected way;
- **routing inside the cores** (`map_and_route`): each core routed on its own
  coupling graph with SABRE, giving a circuit on physical qubits with SWAPs
  and teleportations, optionally through fixed communication qubits (ports);
- **timing**: how long the circuit takes when links make EPR pairs a few at
  a time, and `"hqa-latency"`, an HQA that plans for it.

```python
from qevedo.compiler.multicore import Timing, map_and_route, map_to_cores

mapping = map_to_cores(open("circuit.qasm").read(), cores=4, qubits_per_core=10)
print(mapping.communications)            # teleportations

result = map_and_route(open("circuit.qasm").read(), cores=4, qubits_per_core=9,
                       chip="grid:3x3", topology="line", ports="central",
                       method="hqa-latency", timing=Timing(link_capacity=2))
print(result.teleportations, result.swaps, result.schedule().duration)
```

It reproduces Escofet et al.'s comparison (HQA needs 1.24× fewer
teleportations than FGP-rOEE), and on TeleSABRE's own device and circuits,
HQA with ports uses 0.26× its EPR pairs and 0.51× its SWAPs. See
[docs/MULTICORE.md](https://github.com/qevedo/qcc/blob/master/docs/MULTICORE.md)
for the model, the methods and every benchmark.

## Tests

```bash
pytest
```
