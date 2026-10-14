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
* ``"fgp-roee"``: Fine-Grained Partitioning with relaxed Overall Extreme
  Exchange (Baker et al., 2020), the earlier state of the art.
* ``"naive"``: move one qubit of each split gate to the other's core.

Every method starts from the same initial placement (an OEE partition of how
often each pair of qubits interacts), and the cores are all-to-all inside and
between each other, as in the papers.
"""

from __future__ import annotations

import time
from collections.abc import Sequence
from dataclasses import dataclass

from qevedo.compiler.ir import Circuit, Qubit

__all__ = ["METHODS", "CoreMapping", "map_to_cores", "random_circuit_bounds"]

METHODS = ("hqa", "hqa-no-lookahead", "fgp-roee", "naive")


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
    #: Teleportations: how many times a qubit is in a different core than in
    #: the step before.
    communications: int
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
) -> CoreMapping:
    """Map ``circuit`` onto ``cores`` cores of ``qubits_per_core`` qubits each.

    ``circuit`` is a :class:`~qevedo.compiler.ir.Circuit`, OpenQASM 2 or 3
    source, or the circuit's two-qubit gates as pairs of qubit indices in
    program order (then ``num_qubits`` gives the number of qubits, by default
    one more than the largest index). Gates on three or more qubits are
    decomposed first. Raises ``ValueError`` when the qubits do not fit or the
    method is unknown.
    """
    from qevedo.compiler import _core

    if method not in METHODS:
        raise ValueError(f"unknown mapping method {method!r}; choose one of {', '.join(METHODS)}")
    qubits, pairs = _two_qubit_gates(circuit, num_qubits)
    start = time.perf_counter()
    initial, assignments, timeslices, communications = _core.map_to_cores(
        pairs, len(qubits), cores, qubits_per_core, method
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
        communications=communications,
        seconds=seconds,
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
