"""Mapping circuits onto multi-core (modular) quantum computers.

A multi-core quantum computer is several small processors ("cores") linked by
quantum channels. A two-qubit gate needs both of its qubits in the same core,
and moving a qubit to another core (by teleportation) is slow and noisy. A
mapper decides, before every timeslice of the circuit, which core each qubit
is in, with as few moves (*non-local communications*) as possible.

    from qevedo.compiler.multicore import map_to_cores

    mapping = map_to_cores(open("qft_40.qasm").read(), cores=4, qubits_per_core=10)
    print(mapping.communications, "teleportations")

Methods (all in the Rust core):

* ``"hqa"``: Hungarian Qubit Assignment, from Escofet, Ovide, Almudéver,
  Alarcón and Abadal (arXiv:2309.12182; ACM TQC 2025, arXiv:2403.17205).
  Before each timeslice, the gates split across cores are assigned to cores
  with the Hungarian algorithm, using a cost that counts teleportations and
  favours cores where the qubits will be needed next. The default.
* ``"hqa-latency"``: HQA that also weighs how much later each placement
  lets its gate start, favouring moves of qubits that have been idle: a few
  more teleportations, a shorter circuit (not in the papers).
* ``"hqa-no-lookahead"``: HQA without the look-ahead term.
* ``"hqa-uniform-cost"``: HQA counting every move as one teleportation, as in
  the papers, even when the cores it crosses are not linked directly.
* ``"fgp-roee"``: Fine-Grained Partitioning with relaxed Overall Extreme
  Exchange (Baker et al., 2020), the earlier state of the art.
* ``"naive"``: move one qubit of each split gate to the other's core.

Every method starts from the same initial placement (an OEE partition of how
often each pair of qubits interacts). Cores are all-to-all inside. Between
cores, the links are all-to-all as in the papers, or a line, a ring, a grid or
any connected set of links (``topology``); moving a qubit between cores that
are not linked takes one teleportation per link on the shortest path, and HQA
counts those links in its cost.

When the qubits inside each core are not all coupled to each other (a line,
a grid...), :func:`map_and_route` also routes each core with SABRE and
returns a circuit on physical qubits, with SWAPs inside cores and
teleportations between them:

    from qevedo.compiler.multicore import map_and_route

    result = map_and_route(source, cores=4, qubits_per_core=9, chip="grid:3x3")
    print(result.teleportations, "teleportations,", result.swaps, "SWAPs")
"""

from __future__ import annotations

import time
from collections.abc import Sequence
from dataclasses import astuple, dataclass, replace

from qevedo.compiler.ir import Circuit, Clbit, Instruction, Qubit

__all__ = [
    "METHODS",
    "CoreCircuit",
    "CoreMapping",
    "CoreSchedule",
    "Timing",
    "map_and_route",
    "map_to_cores",
    "random_circuit_bounds",
    "topology_links",
]

METHODS = ("hqa", "hqa-latency", "hqa-no-lookahead", "hqa-uniform-cost", "fgp-roee", "naive")

@dataclass(frozen=True)
class Timing:
    """Durations and link resources, in any one time unit.

    The defaults are illustrative, in units of a two-qubit gate: set them for
    the hardware at hand. Each link between cores has ``link_capacity``
    channels; a channel holds one EPR pair and makes a new one in
    ``epr_pair`` once its pair is used. A teleportation uses one pair per
    link it crosses and then takes ``teleport`` per link.
    """

    one_qubit: float = 0.1
    two_qubit: float = 1.0
    swap: float = 3.0
    measure: float = 5.0
    epr_pair: float = 10.0
    teleport: float = 5.0
    link_capacity: int = 1


#: A ``topology``: a name, or the links between cores as pairs of core indices.
Topology = str | Sequence[tuple[int, int]]


@dataclass
class CoreMapping:
    """Where every qubit is, before the circuit and during each timeslice."""

    method: str
    cores: int
    qubits_per_core: int
    #: The circuit's qubits, in the order their indices refer to.
    qubits: list[Qubit]
    #: ``initial[q]``: the core of qubit ``q`` before the circuit starts.
    initial: list[int]
    #: ``assignments[t][q]``: the core of qubit ``q`` during timeslice ``t``.
    assignments: list[list[int]]
    #: ``timeslices[t]``: the two-qubit gates (qubit index pairs) of timeslice ``t``.
    timeslices: list[list[tuple[int, int]]]
    #: The links between cores, each once with the smaller core first.
    links: list[tuple[int, int]]
    #: Teleportations: every time a qubit is in a different core than in the
    #: step before, one per link on the way.
    communications: int
    #: How many times a qubit is in a different core than in the step before
    #: (equal to ``communications`` when all cores are linked to all).
    moved_qubits: int
    #: Time the mapper took, in seconds.
    seconds: float

    def moves(self) -> list[list[tuple[int, int, int]]]:
        """For each timeslice, the moves before it, as ``(qubit, from_core, to_core)``."""
        out = []
        previous = self.initial
        for assignment in self.assignments:
            out.append([(q, a, b) for q, (a, b) in enumerate(zip(previous, assignment)) if a != b])
            previous = assignment
        return out


def map_to_cores(
    circuit: Circuit | str | Sequence[tuple[int, int]],
    cores: int,
    qubits_per_core: int,
    method: str = "hqa",
    num_qubits: int | None = None,
    topology: Topology = "all-to-all",
    timing: Timing | None = None,
) -> CoreMapping:
    """Map ``circuit`` onto ``cores`` cores of ``qubits_per_core`` qubits each.

    ``circuit`` is a :class:`~qevedo.compiler.ir.Circuit`, OpenQASM 2 or 3
    source, or the circuit's two-qubit gates as pairs of qubit indices in
    program order (then ``num_qubits`` gives the number of qubits, by default
    one more than the largest index). Gates on three or more qubits are
    decomposed first.

    ``topology`` links the cores: ``"all-to-all"`` (the default),
    ``"line"``, ``"ring"``, ``"grid"`` (the squarest grid of ``cores``
    cores), ``"grid:RxC"``, or a list of links as pairs of core indices.

    ``"hqa-latency"`` keeps whichever of its latency weights gives the
    mapping that takes least time under ``timing`` (by default
    :class:`Timing`); the other methods ignore it.

    Raises ``ValueError`` when the qubits do not fit, the method is unknown or
    the links do not connect every core.
    """
    from qevedo.compiler import _core

    if method not in METHODS:
        raise ValueError(f"unknown mapping method {method!r}; choose one of {', '.join(METHODS)}")
    qubits, pairs = _two_qubit_gates(circuit, num_qubits)
    links = topology_links(topology, cores)
    start = time.perf_counter()
    initial, assignments, timeslices, communications, moved = _core.map_to_cores(
        pairs, len(qubits), cores, qubits_per_core, method, links, astuple(timing or Timing())
    )
    seconds = time.perf_counter() - start
    return CoreMapping(
        method=method,
        cores=cores,
        qubits_per_core=qubits_per_core,
        qubits=qubits,
        initial=initial,
        assignments=assignments,
        timeslices=[[tuple(g) for g in s] for s in timeslices],
        links=sorted({(min(a, b), max(a, b)) for a, b in links}),
        communications=communications,
        moved_qubits=moved,
        seconds=seconds,
    )


@dataclass
class CoreCircuit:
    """A circuit mapped onto cores and routed inside them, ready to run."""

    #: Which core each qubit is in, step by step.
    mapping: CoreMapping
    #: The circuit on physical qubits: register ``q`` of ``cores × qubits_per_core``
    #: qubits, where ``q[c · qubits_per_core + s]`` is slot ``s`` of core ``c``.
    #: SWAPs inside cores are ``swap``; moves between cores are ``teleport``
    #: (source, destination). Barriers are dropped.
    circuit: Circuit
    #: How the qubits of each core are coupled, as pairs of slots.
    chip: list[tuple[int, int]]
    #: ``initial_layout[q]``: the physical qubit of circuit qubit ``q`` at the start.
    initial_layout: list[int]
    #: ``final_layout[q]``: where circuit qubit ``q`` ends up.
    final_layout: list[int]
    #: SWAPs inserted inside cores.
    swaps: int
    #: Teleportations, one per link crossed (``mapping.communications``).
    teleportations: int
    #: Time mapping and routing took, in seconds.
    seconds: float
    #: The timing the mapping planned for, and :meth:`schedule`'s default.
    timing: Timing = Timing()

    def schedule(self, timing: Timing | None = None, **changes) -> CoreSchedule:
        """How long the circuit takes when links make EPR pairs at a limited rate.

        Under ``timing`` (by default the one the circuit was mapped for),
        with any field changed by keyword, e.g. ``schedule(link_capacity=2)``
        (see :class:`Timing`). Operations start as soon as their qubits (and,
        for teleportations, an EPR pair) are ready.
        """
        from qevedo.compiler import _core

        raw = [
            (inst.name, [q.index for q in inst.qubits], list(inst.params), [])
            for inst in self.circuit.instructions
        ]
        start, duration, epr_wait, link_hops = _core.schedule(
            raw,
            self.mapping.cores,
            self.mapping.qubits_per_core,
            self.mapping.links,
            astuple(replace(timing or self.timing, **changes)),
        )
        return CoreSchedule(
            start=start,
            duration=duration,
            epr_wait=epr_wait,
            link_hops={tuple(link): hops for link, hops in link_hops},
        )


@dataclass
class CoreSchedule:
    """When each instruction of a :class:`CoreCircuit` runs."""

    #: ``start[i]``: when ``circuit.instructions[i]`` starts.
    start: list[float]
    #: When the last instruction ends.
    duration: float
    #: Total time teleportations waited for an EPR pair, beyond waiting for
    #: their qubit.
    epr_wait: float
    #: How many teleportation hops crossed each link.
    link_hops: dict[tuple[int, int], int]


def map_and_route(
    circuit: Circuit | str | Sequence[tuple[int, int]],
    cores: int,
    qubits_per_core: int,
    chip: Topology = "line",
    method: str = "hqa",
    topology: Topology = "all-to-all",
    num_qubits: int | None = None,
    timing: Timing | None = None,
    trials: int = 8,
    seed: int = 1,
    partner_slots: bool = True,
    look_ahead: int = 4,
    ports: str | dict[tuple[int, int], tuple[int, int]] | None = None,
) -> CoreCircuit:
    """Map ``circuit`` onto cores, as :func:`map_to_cores`, and route it inside them.

    ``circuit`` is a :class:`~qevedo.compiler.ir.Circuit`, OpenQASM source, or
    two-qubit gates as pairs of qubit indices (taken as CX gates on
    ``num_qubits`` qubits, as in :func:`map_to_cores`).

    ``chip`` is how the ``qubits_per_core`` qubits of every core are coupled:
    ``"line"`` (the default), ``"ring"``, ``"grid"``, ``"grid:RxC"``,
    ``"all-to-all"`` or a list of pairs of slots. ``topology`` links the
    cores, as in :func:`map_to_cores`. Between the steps where qubits change
    core, each core's part of the circuit is routed with SABRE (``trials``
    random starts for the first layout); a qubit arriving in a core takes the
    free slot closest to its next partner. ``timing`` is what
    ``"hqa-latency"`` plans for, and what :meth:`CoreCircuit.schedule` uses
    by default.

    ``ports`` fixes where links meet the cores, as on hardware: a qubit must
    then be routed onto its core's port of a link to teleport, and lands on
    the other core's port; the mapper keeps one slot free in every core.
    ``"central"`` puts the ports on the most central slots of each chip; a
    dict maps each link ``(a, b)`` to the slots of its port in cores ``a`` and
    ``b``. By default (``None``), a qubit teleports from any slot to any free
    one.

    For ablation experiments: ``partner_slots=False`` puts an arriving qubit
    in the first free slot instead of next to its next partner, and
    ``look_ahead`` is how many timeslices past each epoch SABRE looks at
    (0: none).

    Raises ``ValueError`` as :func:`map_to_cores` does, and when the chip's
    couplings do not connect its qubits.
    """
    from qevedo.compiler import _core

    if method not in METHODS:
        raise ValueError(f"unknown mapping method {method!r}; choose one of {', '.join(METHODS)}")
    if isinstance(circuit, str):
        from qevedo.compiler.frontend.qasm import parse_qasm

        circuit = parse_qasm(source=circuit)
    if not isinstance(circuit, Circuit):
        qubits, pairs = _two_qubit_gates(circuit, num_qubits)
        circuit = Circuit(
            num_qubits=len(qubits),
            qreg_sizes={"q": len(qubits)},
            instructions=[Instruction("cx", (qubits[a], qubits[b])) for a, b in pairs],
        )
    if any(len(inst.qubits) > 2 and inst.name != "barrier" for inst in circuit.instructions):
        from qevedo.compiler.passes.decompose import DecomposeToNative

        circuit = DecomposeToNative().run(circuit)
    qubits = [Qubit(name, i) for name, size in circuit.qreg_sizes.items() for i in range(size)]
    qubit_index = {q: i for i, q in enumerate(qubits)}
    clbits = [Clbit(name, i) for name, size in circuit.creg_sizes.items() for i in range(size)]
    clbit_index = {c: i for i, c in enumerate(clbits)}
    raw = [
        (
            inst.name,
            [qubit_index[q] for q in inst.qubits],
            list(inst.params),
            [clbit_index[c] for c in inst.clbits],
        )
        for inst in circuit.instructions
    ]
    links = topology_links(topology, cores)
    chip_edges = topology_links(chip, qubits_per_core)
    port_slots = _port_slots(ports, links)
    start = time.perf_counter()
    (
        initial,
        assignments,
        timeslices,
        communications,
        moved,
        physical,
        initial_layout,
        final_layout,
        swaps,
        teleportations,
    ) = _core.map_and_route(
        raw,
        len(qubits),
        cores,
        qubits_per_core,
        method,
        links,
        chip_edges,
        astuple(timing or Timing()),
        trials,
        seed,
        partner_slots,
        look_ahead,
        port_slots,
        ports == "central",
    )
    seconds = time.perf_counter() - start
    mapping = CoreMapping(
        method=method,
        cores=cores,
        qubits_per_core=qubits_per_core,
        qubits=qubits,
        initial=initial,
        assignments=assignments,
        timeslices=[[tuple(g) for g in s] for s in timeslices],
        links=sorted({(min(a, b), max(a, b)) for a, b in links}),
        communications=communications,
        moved_qubits=moved,
        seconds=seconds,
    )
    out = Circuit(
        num_qubits=cores * qubits_per_core,
        num_clbits=circuit.num_clbits,
        qreg_sizes={"q": cores * qubits_per_core},
        creg_sizes=dict(circuit.creg_sizes),
        instructions=[
            Instruction(
                name,
                tuple(Qubit("q", p) for p in qs),
                tuple(params),
                tuple(clbits[c] for c in cs),
            )
            for name, qs, params, cs in physical
        ],
    )
    return CoreCircuit(
        mapping=mapping,
        circuit=out,
        chip=sorted({(min(a, b), max(a, b)) for a, b in chip_edges}),
        initial_layout=initial_layout,
        final_layout=final_layout,
        swaps=swaps,
        teleportations=teleportations,
        seconds=seconds,
        timing=timing or Timing(),
    )


def _port_slots(
    ports: str | dict[tuple[int, int], tuple[int, int]] | None, links: list[tuple[int, int]]
) -> list[tuple[int, int]] | None:
    """Port slots per link, in the core's order of links (sorted, smaller core first)."""
    if ports is None or ports == "central":
        return None
    if isinstance(ports, str):
        raise ValueError(f"unknown ports {ports!r}: use 'central' or a dict of links to slots")
    out = []
    for a, b in sorted({(min(a, b), max(a, b)) for a, b in links}):
        if (a, b) in ports:
            out.append(tuple(ports[(a, b)]))
        elif (b, a) in ports:
            sb, sa = ports[(b, a)]
            out.append((sa, sb))
        else:
            raise ValueError(f"no port for the link ({a}, {b})")
    return out


def topology_links(topology: Topology, cores: int) -> list[tuple[int, int]]:
    """The links between ``cores`` cores for a ``topology`` name or list of links.

    The same names describe how the qubits inside a core are coupled.
    """
    if not isinstance(topology, str):
        return [(int(a), int(b)) for a, b in topology]
    name, _, shape = topology.partition(":")
    if name == "all-to-all":
        return [(a, b) for a in range(cores) for b in range(a + 1, cores)]
    if name == "line":
        return [(c - 1, c) for c in range(1, cores)]
    if name == "ring":
        return [(c - 1, c) for c in range(1, cores)] + ([(cores - 1, 0)] if cores > 2 else [])
    if name == "grid":
        if shape:
            rows, columns = (int(x) for x in shape.lower().split("x"))
        else:
            rows = max(r for r in range(1, int(cores**0.5) + 1) if cores % r == 0)
            columns = cores // rows
        if rows * columns != cores:
            raise ValueError(f"a {rows}x{columns} grid has {rows * columns} cores, not {cores}")
        links = []
        for r in range(rows):
            for c in range(columns):
                core = r * columns + c
                if c + 1 < columns:
                    links.append((core, core + 1))
                if r + 1 < rows:
                    links.append((core, core + columns))
        return links
    raise ValueError(
        f"unknown topology {topology!r}: use all-to-all, line, ring, grid, grid:RxC or a list of links"
    )


def random_circuit_bounds(cores: int, qubits: int, two_qubit_gates: int) -> tuple[float, float]:
    """Expected lower and upper bounds on communications for a random circuit.

    From Escofet et al. (ACM TQC 2025, Theorems 3.2 and 3.3): with gates on
    random pairs of ``qubits`` qubits spread over ``cores`` cores, a mapper
    that leaves qubits where they are needs about ``(N − 1) G q / (N (q − 1))``
    teleportations for ``G`` two-qubit gates, and the naive strategy about
    twice that. Mappers that move qubits ahead of time (HQA) can beat the
    lower one.
    """
    from qevedo.compiler import _core

    return _core.random_circuit_bounds(cores, qubits, two_qubit_gates)


def _two_qubit_gates(
    circuit: Circuit | str | Sequence[tuple[int, int]], num_qubits: int | None
) -> tuple[list[Qubit], list[tuple[int, int]]]:
    """The circuit's qubits and its two-qubit gates as index pairs."""
    if isinstance(circuit, str):
        from qevedo.compiler.frontend.qasm import parse_qasm

        circuit = parse_qasm(source=circuit)
    if not isinstance(circuit, Circuit):
        pairs = [(int(a), int(b)) for a, b in circuit]
        n = num_qubits if num_qubits is not None else 1 + max((max(p) for p in pairs), default=-1)
        return [Qubit("q", i) for i in range(n)], pairs
    if any(len(inst.qubits) > 2 and inst.name not in ("barrier",) for inst in circuit.instructions):
        from qevedo.compiler.passes.decompose import DecomposeToNative

        circuit = DecomposeToNative().run(circuit)
    qubits = [Qubit(name, i) for name, size in circuit.qreg_sizes.items() for i in range(size)]
    index = {q: i for i, q in enumerate(qubits)}
    pairs = [
        (index[inst.qubits[0]], index[inst.qubits[1]])
        for inst in circuit.instructions
        if len(inst.qubits) == 2 and inst.name not in ("barrier", "measure", "reset")
    ]
    return qubits, pairs
