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
from qevedo.compiler.multicore import Timing, map_and_route, map_to_cores

mapping = map_to_cores(open("circuit.qasm").read(), cores=4, qubits_per_core=10)
print(mapping.communications)    # teleportations needed
print(mapping.assignments[0])    # core of each qubit during the first timeslice
print(mapping.moves()[3])        # (qubit, from core, to core) before timeslice 3

# Cores in a line, a ring, a grid, or linked by any connected set of links:
mapping = map_to_cores(circuit, cores=8, qubits_per_core=10, topology="grid:2x4")
mapping = map_to_cores(circuit, cores=4, qubits_per_core=10, topology=[(0, 1), (1, 2), (1, 3)])

# Qubits inside each core coupled in a 3 × 3 grid: route each core with SABRE
# and get a circuit on physical qubits, with SWAPs and teleportations.
result = map_and_route(circuit, cores=4, qubits_per_core=9, chip="grid:3x3")
print(result.teleportations, result.swaps)
print(result.circuit.instructions[:5])   # cx, swap and teleport on q[0..35]

# How long it takes when each link holds two EPR pairs at a time.
print(result.schedule(link_capacity=2).duration)

# Map for time rather than teleportations, for links with two channels.
fast = map_and_route(circuit, cores=4, qubits_per_core=9, chip="grid:3x3",
                     method="hqa-latency", timing=Timing(link_capacity=2))
```

`map_to_cores` also takes a parsed `Circuit` or a list of two-qubit gates as
qubit pairs, and a `method`:

| method | what it does |
|---|---|
| `"hqa"` (default) | Hungarian Qubit Assignment with look-ahead attraction |
| `"hqa-latency"` | HQA that also weighs how much each move delays its gate, planning for a `timing` (not in the papers; see below) |
| `"hqa-no-lookahead"` | HQA without the attraction term |
| `"hqa-uniform-cost"` | HQA counting every move as one teleportation (the papers' cost), even across unlinked cores |
| `"fgp-roee"` | Fine-Grained Partitioning with relaxed Overall Extreme Exchange |
| `"naive"` | move one qubit of each split gate to the other's core |

From Rust, `qcc::multicore::{Architecture, Timeslices, map, Method}` (and
`map_and_route` with a `qcc::routing::CouplingGraph` for the chip, and
`map_timed` with a `Timing` for `HqaLatency`; `hqa` with `HqaOptions` sets
the latency weight and the teleportation time directly), with
`Architecture::new` (all to all), `line`, `ring`, `grid` or `with_links` (see
`core/examples/multicore.rs`, run with
`cargo run --release --example multicore -- 100 10`).

## The model

As in the papers:

- `cores` cores of `qubits_per_core` physical qubits each;
- each core is all-to-all inside, and the cores are all-to-all between each
  other, so moving a qubit's state to any other core is one communication
  (beyond the papers, the cores can instead be linked in a line, a ring, a
  grid or any connected way; then a move takes one communication per link on
  the shortest path);
- the circuit is cut into **timeslices**, sets of two-qubit gates that run in
  parallel (each qubit in at most one gate per slice), built as soon as
  possible; single-qubit gates do not constrain the mapping;
- the cost is the number of times a qubit is in a different core than in the
  step before (times the links crossed, on sparse links). The initial
  placement is free.

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

On sparse links, the papers note the cost "can be adapted"; here both terms
count teleportations: the links qA and qB cross to reach c, plus, for the
look-ahead, the links each would cross from c to meet its future partners,
weighted by `2^−(m−t)`. With every core linked to every other this differs
from the formula above by a constant per gate, so the Hungarian algorithm
makes the same choices.

## Routing inside the cores

The mappers, like the papers, take every core as all-to-all inside. When the
qubits of a core are coupled in a line, a grid or any other graph (the
*chip*), `map_and_route` also produces the circuit that runs:

1. The timeslices are grouped into **epochs**, runs in which no qubit changes
   core.
2. In each epoch, each core runs its part of the circuit, single-qubit gates
   and measurements included, routed on its chip with **SABRE** (Li, Ding and
   Xie, ASPLOS 2019), which inserts SWAPs. Each core starts where the
   previous epoch left its qubits; the first epoch's layout comes from
   SABRE's forward–backward layout search. SABRE also looks at the gates
   the core's qubits have together in the next four timeslices, so an
   epoch ends with them ready.
3. Between epochs, the qubits that change core are **teleported**. A qubit
   arriving in a core takes the free slot closest to the qubit it interacts
   with next, or, if that one has not arrived yet, the free slot with most
   free neighbours, so the partner can land next to it.

The output is a circuit on `cores × qubits_per_core` physical qubits, where
`q[c · qubits_per_core + s]` is slot `s` of core `c`, with `swap` inside
cores and `teleport source, destination` between them. The teleportations
before a timeslice happen at once (as with communication qubits in each
core): a destination is empty or is being left by another teleportation.
They are listed so that each destination is already free whenever possible.
Tests simulate the routed circuits and check that they compute the same state
as the original.

The SABRE router is a module of its own (`core/src/routing/`) and routes any
circuit on any connected chip; on 20-qubit lines, 5 × 5 grids and 19-qubit
heavy-hex chips it inserts about as many SWAPs as Qiskit's `SabreLayout`
(2774 against 2806 over QFT, random, Quantum Volume and Cuccaro circuits;
`benchmarks/routing/vs_qiskit.py`, results in
`benchmarks/routing/results/sabre.csv`). It is slower than Qiskit's, which
is also written in Rust: 5–110 ms per circuit against 2–21 ms.

## Time on links of limited capacity

Counting teleportations treats every one alike, but a link makes EPR pairs
at a limited rate, a few at a time: when many qubits cross the same link
together, they wait. `CoreCircuit.schedule` (Rust: `multicore::schedule`)
schedules the routed circuit as soon as possible under a simple model:

- each operation takes a fixed time by kind (single- and two-qubit gates,
  SWAPs, measurements);
- each link has `link_capacity` channels; a channel holds one EPR pair and
  makes a new one in `epr_pair` once its pair is used;
- a teleportation uses one pair per link it crosses, one after the other,
  and then takes `teleport` per link; the teleportations before a timeslice
  run at once, as long as there are pairs.

The defaults (in units of a two-qubit gate: SWAP 3, measurement 5, EPR pair
10, the rest of a teleportation 5, one channel) are illustrative; set them
for the hardware at hand. The result has the total duration, each
instruction's start, the time teleportations waited for pairs and the hops
per link.

## Code

Rust, in `core/src/multicore/`, one concept per file:

| file | what it holds |
|---|---|
| `mod.rs` | overview, `Method`, `initial_placement`, `map` |
| `architecture.rs` | the hardware model: cores and the links between them |
| `timeslices.rs` | circuits as timeslices of two-qubit gates |
| `assignment.rs` | assignments, mappings, validity, communication count |
| `lookahead.rs` | how future interactions attract qubits to cores |
| `hungarian.rs` | the Hungarian algorithm (minimum-cost assignment) |
| `oee.rs` | Overall Extreme Exchange graph partitioning |
| `hqa.rs` | Hungarian Qubit Assignment |
| `baselines.rs` | the naive and FGP-rOEE mappers |
| `physical.rs` | routing inside cores and teleportations: a circuit that runs |
| `schedule.rs` | timing: how long that circuit takes on links of limited capacity |
| `bounds.rs` | theoretical bounds for random circuits |

and in `core/src/routing/`, `coupling.rs` (chips: coupling graphs and their
distances) and `sabre.rs` (the SABRE router and layout search).

Python: `python/qevedo/compiler/multicore.py` (API) and the
`map_to_cores`, `map_and_route` and `schedule` bindings in `python/src/lib.rs`. Tests: the `tests` module of
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

- **Sparse links** (not in the papers): the cost counts link hops now and in
  the look-ahead, as above. A small cost per moved qubit (0.01) breaks ties
  between moving one qubit and moving both to a core in between, which cost
  as many hops but disturb two placements; it was chosen on these benchmarks,
  but any value from 0.001 to 0.05 gives totals within 1%. The initial OEE
  groups go to cores so that groups that interact most are close.

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
| QFT | 1.84× | 1.23× | 20.9× |
| Draper adder | 1.44× | 3.74× | 10.8× |
| Cuccaro adder | 0.73× | 4.30× | 4.9× |
| Random | 1.10× | 1.87× | 6.5× |
| Quantum Volume | 1.16× | 1.39× | 9.2× |
| **all points** | **1.24× (mean)** | **3.45× (mean)** | **7.8× (median)** |

Medians per benchmark over the virtual, weak and strong scaling
configurations; above 1 means HQA needs fewer communications or less time.

- HQA needs **1.24× fewer communications** than FGP-rOEE on average, in line
  with the 1.28× of the original HQA paper (the extended paper reports 1.56×
  over a different set of configurations), and is a median **8× faster**.
- The look-ahead attraction matters as the paper says: on cores of 16
  qubits, as many as needed, HQA without it needs 4.04× as many
  communications on Cuccaro adders of 20 to 200 qubits (1.32–8.86×) and
  1.18× on random circuits, where the 2025 paper reports 2.27× and 1.09×
  (`results/attraction.csv`).
- On the Cuccaro adder, whose gates mostly join neighbouring qubits,
  FGP-rOEE needs fewer communications than HQA here (0.73×).
- On random circuits HQA goes below the papers' "lower bound", which assumes
  qubits only move when a gate needs them; HQA moves them ahead of time.

### Sparse links

An experiment the papers do not have: 100-qubit circuits on 10 cores of 10
qubits, linked all to all, in a ring, in a 2 × 5 grid or in a line.
Communications count link hops. Mean ratios over the five benchmarks:

| links | papers' cost / HQA | FGP-rOEE / HQA | naive / HQA |
|---|---|---|---|
| all to all | 1.00× | 1.31× | 2.34× |
| ring | 1.06× | 1.91× | 3.84× |
| 2 × 5 grid | 1.04× | 1.60× | 3.39× |
| line | 1.10× | 1.76× | 3.17× |

- Counting link hops in HQA's cost saves 4–10% over the papers' cost, most on
  the line, where distances are largest.
- HQA's lead over FGP-rOEE grows on sparse links: 1.6–1.9× fewer
  communications, against 1.3× with all cores linked to all.

### Routing inside the cores

The same 100-qubit circuits on 10 cores of 10 qubits (cores linked all to
all), now with the qubits of each core coupled in a 2 × 5 grid or a line
(`results/routing.csv`). Totals over the five benchmarks:

| chip | method | teleportations | SWAPs |
|---|---|---|---|
| 2 × 5 grid | HQA | **10193** | 13628 |
| | FGP-rOEE | 13023 | 12411 |
| | naive | 17790 | 10106 |
| line | HQA | **10193** | 26943 |
| | FGP-rOEE | 13023 | 25028 |
| | naive | 17790 | 23480 |

- Fewer teleportations can cost SWAPs: the naive mapper moves a qubit only
  to meet its partner, so it lands next to it, while HQA also moves qubits
  ahead of time, before their partners are there.
- Still, HQA has the lowest total wherever one teleportation costs more than
  1.4 SWAPs (the largest break-even over every benchmark and both chips),
  and a teleportation, which consumes an EPR pair shared across a link and
  needs a Bell measurement, costs far more than a local SWAP on any hardware
  proposed. The exception is again the Cuccaro adder, where FGP-rOEE needs
  fewer teleportations (114 against 156), and on a line also fewer SWAPs.
- Choosing the slot of an arriving qubit by its next partner, rather than
  the first free slot, saves 18% of the SWAPs over all these runs; looking
  at the next timeslices saves 10% more (any horizon from 1 to 64 slices
  gives totals within 2%). `results/ablation.csv` has each variant
  (`map_and_route(..., partner_slots=False, look_ahead=0)`).

### Duration with links of limited capacity

The same circuits routed on 2 × 5-grid cores, cores linked all to all or in
a line, scheduled with the default timing and 1, 2, 4 or unlimited channels
per link (`results/latency.csv`). Mean, over the five benchmarks, of each
method's duration over HQA's (above 1: HQA finishes first):

| links | channels | HQA-latency / HQA | FGP-rOEE / HQA | naive / HQA |
|---|---|---|---|---|
| all to all | 1 | **0.94×** | 1.21× | 1.64× |
| | 2 | **0.93×** | 1.00× | 1.27× |
| | 4 | **0.93×** | 0.97× | 1.06× |
| | unlimited | **0.92×** | 0.96× | 1.03× |
| line | 1 | **0.91×** | 1.49× | 2.25× |
| | 2 | **0.91×** | 1.38× | 2.08× |
| | 4 | **0.90×** | 1.27× | 1.82× |
| | unlimited | **0.88×** | 1.20× | 1.51× |

- When pairs are scarce, fewer teleportations means a shorter circuit:
  with one channel per link HQA finishes 1.21× sooner than FGP-rOEE with
  cores linked all to all, and 1.5× sooner on a line, where every crossing
  shares the few links (except, again, the Cuccaro adder).
- When links have pairs to spare, teleportations run in parallel and the
  count matters less: with all cores linked to all and four channels or
  more, FGP-rOEE and even the naive mapper finish about as soon as HQA, or
  slightly sooner. HQA minimizes teleportations, not the length of the
  critical path.
- **HQA-latency** (`"hqa-latency"`, not in the papers) closes that gap. HQA
  keeps, for every qubit, when it would be free if everything ran as soon as
  possible, and adds to the cost of putting a gate in a core how much later
  the gate could start than if nothing moved, in teleportations, times a
  weight. Between moving one qubit or the other, it moves the one that has
  been idle, whose teleportation overlaps other gates. No single weight
  suits every circuit: weighing delay shortens circuits whose links have
  pairs to spare, but its extra teleportations lengthen those that keep the
  links busy (with a fixed weight of 0.5, random circuits and Quantum Volume
  on a line with few channels took 8–18% longer than with HQA). So
  HQA-latency maps with the weights 0, 0.25, 0.5 and 1 and keeps the one
  that runs fastest under the `Timing` it is given: `map_and_route` routes
  all four and schedules them; `map_to_cores`, which does not route, uses an
  estimate without SWAPs (gates and teleportations as soon as possible,
  sharing the links' channels; `mapping_duration`). Weight 0 is plain HQA,
  so through `map_and_route` it is never slower than HQA: 6–12% shorter on
  average, up to 25% (QFT on a line), the shortest of all four methods in
  31 of the 40 runs. It uses about as many teleportations as HQA when
  channels are few, and up to 11% more when there are pairs to spare. The
  exceptions are Cuccaro (FGP-rOEE, as before), random circuits with two channels (FGP-rOEE, by 1%) and, with all
  cores linked and four channels or more, QFT, Draper and Quantum Volume,
  where the naive mapper is 2–11% shorter. Mapping and routing take about
  twice as long as with HQA (0.08 s against 0.03 s for 100-qubit Quantum
  Volume).
- Making the delay term also count how busy each link is (when its next
  EPR pair is ready) gave longer circuits, not shorter, so it was left out.

#### Faster and slower links

The default timing is illustrative, so `results/timing.csv` repeats the
experiment with links faster (EPR pair and teleportation 2) and slower (50
and 20). Mean duration over HQA's:

| links | timing | channels | HQA-latency | FGP-rOEE | naive |
|---|---|---|---|---|---|
| all to all | fast | 1 | 0.94× | 0.99× | 1.01× |
| | fast | unlimited | 0.94× | 0.97× | 0.96× |
| | slow | 1 | 0.94× | 1.33× | 2.10× |
| | slow | unlimited | 0.88× | 0.92× | 1.20× |
| line | fast | 1 | 0.98× | 1.28× | 1.75× |
| | fast | unlimited | 0.96× | 1.15× | 1.27× |
| | slow | 1 | 0.95× | 1.56× | 2.56× |
| | slow | unlimited | 0.82× | 1.19× | 1.82× |

HQA-latency is again never slower than HQA, 2–18% faster on average; how
much fewer teleportations matter depends on the links: with fast links all
to all, FGP-rOEE and the naive mapper are about as fast as HQA.

### Communication through ports

By default a qubit teleports from any slot of its core to any free slot of
another, as if every qubit could reach every link. With ports
(`map_and_route(..., ports="central")`, or a dict giving, for each link
`(a, b)`, the slots of its port in cores `a` and `b`; `Ports` in Rust), each
link meets each core at one physical qubit, as on hardware:

- to cross a link, a qubit is routed with SWAPs onto its core's port and
  lands on the other core's port, which is emptied first by moving the
  nearest free slot there; a move between unlinked cores crosses the links
  one by one;
- the mapper keeps one slot free in every core (as TeleSABRE does), so a
  core can take a qubit in before another leaves;
- the moves before a timeslice run one after another, each when its next
  core, on any shortest path, has room. When every next core is full (two
  qubits crossing between full cores), two qubits exchange places across a
  link, a teleportation each way at once; if they are not both moving that
  way, one moves back later, costing teleportations the mapping did not
  count (none in a test of 60 random machines filled to capacity).

The code is in `core/src/multicore/ports.rs`.

### Against TeleSABRE

TeleSABRE (Russo et al., IEEE QCE 2025) routes onto multi-core chips with
SWAPs, teleportations and remote gates (telegates) through communication
qubits at fixed places. `benchmarks/telesabre/compare.py` runs its reference
implementation (default configuration) on its device H, 6 cores of 4 × 4
qubits in a 2 × 3 grid of cores, and maps the same circuits onto the same
machine with `map_and_route`, with ports at device H's communication qubits.
EPR pairs / SWAPs (`results/telesabre.csv`; * TeleSABRE's MQT Bench circuits):

| circuit | TeleSABRE | HQA | HQA-latency | FGP-rOEE | HQA, abstract |
|---|---|---|---|---|---|
| Amplitude estimation* | 595 / 2321 | 217 / 1451 | 344 / 1967 | 376 / 1656 | 214 |
| GHZ* | 116 / 421 | 15 / 98 | 16 / 103 | 24 / 136 | 13 |
| Graph state* | 92 / 382 | 11 / 102 | 11 / 102 | 19 / 138 | 30 |
| QFT* | 514 / 2072 | 266 / 1687 | 330 / 1879 | 386 / 1862 | 370 |
| QNN* | 1560 / 6841 | 792 / 5162 | 1073 / 6183 | 1339 / 8317 | 741 |
| Random* | fails | 801 / 4072 | 838 / 4255 | 842 / 4466 | 787 |
| QFT | 1053 / 4163 | 496 / 3085 | 458 / 2509 | 857 / 4896 | 379 |
| Draper adder | 641 / 2861 | 362 / 2367 | 362 / 2367 | 515 / 3356 | 344 |
| Cuccaro adder | 202 / 813 | 12 / 242 | 12 / 242 | 16 / 267 | 38 |
| Random | fails | 1540 / 7603 | 1540 / 7603 | 1594 / 8065 | 1593 |
| Quantum Volume | fails | 2224 / 10235 | 2512 / 11176 | 2576 / 11876 | 2194 |

- TeleSABRE deadlocks on every one of its 100 attempts on the random
  circuits and Quantum Volume. On the other eight, HQA with ports uses 0.26×
  as many EPR pairs and 0.51× as many SWAPs (geometric means), fewer of both
  on every circuit. Mapping and routing take at most 30 ms per circuit
  (parsing aside), against 0.1–43 s for TeleSABRE's completed runs,
  retries included.
- TeleSABRE's paper reports 28% fewer inter-core operations than HQA on a
  similar machine; how HQA's mapping was turned into operations there is
  not described in enough detail to reproduce.
- Two differences remain: TeleSABRE needs its communication qubits free
  when they teleport, while our ports can hold data in between, and our
  exchanges between full cores assume a teleportation each way at once.
- Ports cost many SWAPs: compare with the last column, where communication
  is abstract. HQA-latency minimizes time, not EPR pairs or SWAPs.

#### Time through ports

`results/latency_ports.csv` repeats the latency experiment through central
ports (90-qubit circuits, since the mapper keeps a slot free per core).
Mean duration over HQA's:

| links | channels | HQA-latency | FGP-rOEE | naive |
|---|---|---|---|---|
| all to all | 1 | 1.00× | 1.41× | 2.13× |
| | unlimited | 1.00× | 1.37× | 2.04× |
| line | 1 | 0.93× | 1.66× | 2.80× |
| | unlimited | 0.93× | 1.57× | 2.72× |

Moving qubits to and from the ports makes teleportations dearer, and HQA's
lead grows. The number of channels hardly matters: with one port per link
in each core, the qubits crossing a link already queue at its port.
Through ports, HQA-latency also tries its weights with a teleportation made
dearer by the mean SWAPs to a port and from the other (and keeps the fastest
routed circuit of all eight). That shortens circuits by 7% on a line of
cores, but not at all with every core linked to every other: there the delay
term still misses how crossings queue at the ports.

![Sparse links](../benchmarks/multicore/charts/topology.svg)
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
- G. Li, Y. Ding and Y. Xie, "Tackling the Qubit Mapping Problem for
  NISQ-Era Quantum Devices", ASPLOS (2019),
  [arXiv:1809.02573](https://arxiv.org/abs/1809.02573): SABRE.
- A. Ovide et al., "Mapping quantum algorithms to multi-core quantum computing
  architectures", [arXiv:2303.16125](https://arxiv.org/abs/2303.16125) (2023).
