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
"""

from __future__ import annotations

import time
from collections.abc import Sequence
from dataclasses import dataclass

from qevedo.compiler.ir import Circuit, Qubit

__all__ = ["METHODS", "CoreMapping", "map_to_cores", "random_circuit_bounds", "topology_links"]

METHODS = ("hqa", "hqa-no-lookahead", "hqa-uniform-cost", "fgp-roee", "naive")

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
        pairs, len(qubits), cores, qubits_per_core, method, links
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


def topology_links(topology: Topology, cores: int) -> list[tuple[int, int]]:
    """The links between ``cores`` cores for a ``topology`` name or list of links."""
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
