# Mapping onto multi-core quantum computers

Quantum processors are hard to grow as single chips, so one route to scale is
a **modular** or **multi-core** machine: several small processors (cores)
connected by quantum links. A two-qubit gate can only run when both qubits are
in the same core, and moving a qubit to another core (by teleportation, using
an EPR pair) is slow and noisy. Before every step of the circuit, a *mapper*
decides which core each qubit is in, with as few moves, or **non-local
communications**, as possible.

`qevedo-compiler` implements the Hungarian Qubit Assignment (HQA) mapper of
Escofet et al. and the mappers it was compared with, in the Rust core, with a
Python API.

## Quick start

```python
from qevedo.compiler.multicore import map_to_cores

mapping = map_to_cores(open("circuit.qasm").read(), cores=4, qubits_per_core=10)
print(mapping.communications)    # teleportations needed
print(mapping.assignments[0])    # core of each qubit during the first timeslice
print(mapping.moves()[3])        # (qubit, from core, to core) before timeslice 3
```

`map_to_cores` also takes a parsed `Circuit` or a list of two-qubit gates as
qubit pairs, and a `method`:

| method | what it does |
|---|---|
| `"hqa"` (default) | Hungarian Qubit Assignment with look-ahead attraction |
| `"hqa-no-lookahead"` | HQA without the attraction term |
| `"fgp-roee"` | Fine-Grained Partitioning with relaxed Overall Extreme Exchange |
| `"naive"` | move one qubit of each split gate to the other's core |

From Rust, `qcc::multicore::{Architecture, Timeslices, map, Method}` (see
`core/examples/multicore.rs`, run with
`cargo run --release --example multicore -- 100 10`).

## The model

As in the papers:

- `cores` cores of `qubits_per_core` physical qubits each;
- each core is all-to-all inside, and the cores are all-to-all between each
  other, so moving a qubit's state to any other core is one communication;
- the circuit is cut into **timeslices**, sets of two-qubit gates that run in
  parallel (each qubit in at most one gate per slice), built as soon as
  possible; single-qubit gates do not constrain the mapping;
- the cost is the number of times a qubit is in a different core than in the
  step before. The initial placement is free.

## How HQA works

Before each timeslice, the gates whose qubits are in different cores
("unfeasible" gates) are collected and their qubits taken out of their cores.
Gates, not qubits, are then assigned to cores: the Hungarian algorithm solves
a gates × cores cost matrix, giving each core at most one gate per round, and
both qubits of a gate go into the free space of its core. Rounds repeat until
every gate is placed. The cost of putting gate (qA, qB) in core c is

```
∞                if c has fewer than two free slots
1 − attr(op, c)  if qA or qB was in c before the slice (one teleportation)
2 − attr(op, c)  otherwise (two teleportations)
```

where `attr(op, c)` is the mean attraction of qA and qB to c: the sum, over
their future interactions, of `2^−(m−t)` for interactions `m − t` slices
ahead with a qubit now in c. Every mapper starts from the same placement: an
OEE partition of how often each pair of qubits interacts.

## Code

Rust, in `core/src/multicore/`, one concept per file:

| file | what it holds |
|---|---|
| `mod.rs` | overview, `Method`, `initial_placement`, `map` |
| `architecture.rs` | the hardware model |
| `timeslices.rs` | circuits as timeslices of two-qubit gates |
| `assignment.rs` | assignments, mappings, validity, communication count |
| `lookahead.rs` | how future interactions attract qubits to cores |
| `hungarian.rs` | the Hungarian algorithm (minimum-cost assignment) |
| `oee.rs` | Overall Extreme Exchange graph partitioning |
| `hqa.rs` | Hungarian Qubit Assignment |
| `baselines.rs` | the naive and FGP-rOEE mappers |
| `bounds.rs` | theoretical bounds for random circuits |

Python: `python/qevedo/compiler/multicore.py` (API) and the
`map_to_cores` binding in `python/src/lib.rs`. Tests: the `tests` module of
each Rust file (including a randomized test of every method on 300 random
architectures and circuits) and `python/tests/test_multicore.py`.

## Choices beyond the papers

The papers leave some details open, and use cores of even size. Where the
code had to decide, it says so in `hqa.rs` and `baselines.rs`:

- **Ties** in the Hungarian algorithm go to the lower-numbered gate and core.
- **Look-ahead** stops after 64 slices, past which `2^−64` terms are below
  rounding.
- **Free space** must come in pairs. As in the papers, two cores with an odd
  number of free slots get an auxiliary gate between an idle qubit of each
  (the one least attracted to its core), only when the free space would
  otherwise not hold the gates.
- **Odd core sizes**: a timeslice can then hold more gates than the cores can
  run at once, so `Timeslices::fitted_to` splits it. If the auxiliary-gate rule
  cannot apply (no idle qubit in the odd cores), an idle qubit moves between
  them, or, as a last resort, every idle qubit is taken out and put back after
  the gates; the naive mapper and FGP-rOEE likewise fall back to moving the
  other qubit, or both, when the first move is impossible.

## Reproducing the papers' comparison

`benchmarks/multicore/` runs the experiments of the 2025 paper on its
benchmarks (QFT, Draper and Cuccaro adders, random circuits of depth twice the
qubits, Quantum Volume; built with Qiskit and transpiled to CX):

```bash
pip install qiskit matplotlib
python benchmarks/multicore/run.py      # about a minute; --quick for a few points
python benchmarks/multicore/plot.py     # charts and the summary table
```

Results (`benchmarks/multicore/results/multicore.csv`, Intel i9-14900K):

| benchmark | FGP-rOEE / HQA communications | naive / HQA | FGP-rOEE / HQA time |
|---|---|---|---|
| QFT | 1.84× | 1.23× | 22.2× |
| Draper adder | 1.44× | 3.74× | 13.2× |
| Cuccaro adder | 0.73× | 4.30× | 4.2× |
| Random | 1.10× | 1.87× | 5.6× |
| Quantum Volume | 1.16× | 1.39× | 9.4× |
| **mean over all points** | **1.24×** | **3.45×** | median 7.1× |

Medians per benchmark over the virtual, weak and strong scaling
configurations; above 1 means HQA needs fewer communications or less time.

- HQA needs **1.24× fewer communications** than FGP-rOEE on average, in line
  with the 1.28× of the original HQA paper (the extended paper reports 1.56×
  over a different set of configurations), and is a median **7× faster**.
- The look-ahead attraction matters as the paper says: on Cuccaro adders with
  16 qubits per core, HQA without it needs 1.7–3.9× more communications
  (the paper: 2.27× on average).
- On the Cuccaro adder, whose gates mostly join neighbouring qubits,
  FGP-rOEE needs fewer communications than HQA here (0.73×).
- On random circuits HQA goes below the papers' "lower bound", which assumes
  qubits only move when a gate needs them; HQA moves them ahead of time.

![Strong scaling](../benchmarks/multicore/charts/strong.svg)
![Weak scaling](../benchmarks/multicore/charts/weak.svg)
![Virtual scaling](../benchmarks/multicore/charts/virtual.svg)

## References

- P. Escofet, A. Ovide, C. G. Almudéver, E. Alarcón and S. Abadal, "Hungarian
  Qubit Assignment for Optimized Mapping of Quantum Circuits on Multi-Core
  Architectures", [arXiv:2309.12182](https://arxiv.org/abs/2309.12182) (2023).
- P. Escofet et al. (with M. Bandic, L. Prielinger, H. van Someren, S. Feld
  and C. G. Almudéver), "Revisiting the Mapping of Quantum Circuits: Entering
  the Multi-Core Era", ACM Transactions on Quantum Computing 6(1) (2025),
  [arXiv:2403.17205](https://arxiv.org/abs/2403.17205).
- J. M. Baker et al., "Time-sliced quantum circuit partitioning for modular
  architectures", Computing Frontiers (2020): FGP-rOEE.
- A. Ovide et al., "Mapping quantum algorithms to multi-core quantum computing
  architectures", [arXiv:2303.16125](https://arxiv.org/abs/2303.16125) (2023).
