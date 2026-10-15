"""Mapping onto multi-core architectures (qevedo.compiler.multicore)."""

import random

import numpy as np
import pytest

from qevedo.compiler.multicore import (
    METHODS,
    Timing,
    map_and_route,
    map_to_cores,
    random_circuit_bounds,
    topology_links,
)
from qevedo.compiler.passes import native
from qevedo.compiler.synthesis.gates import gate_matrix

pytestmark = pytest.mark.skipif(not native.available(), reason="the Rust core is not built")


def random_pairs(qubits, layers, seed):
    rng = random.Random(seed)
    pairs = []
    for _ in range(layers):
        order = list(range(qubits))
        rng.shuffle(order)
        pairs += [(order[2 * k], order[2 * k + 1]) for k in range(qubits // 2)]
    return pairs


def check(mapping):
    """Every gate has both qubits in one core and no core is over capacity."""
    for assignment, gates in zip(mapping.assignments, mapping.timeslices):
        for a, b in gates:
            assert assignment[a] == assignment[b]
        for core in range(mapping.cores):
            assert assignment.count(core) <= mapping.qubits_per_core
    moved = sum(len(step) for step in mapping.moves())
    assert moved == mapping.moved_qubits


@pytest.mark.parametrize("method", METHODS)
def test_every_method_gives_a_valid_mapping(method):
    pairs = random_pairs(40, 30, seed=1)
    mapping = map_to_cores(pairs, cores=4, qubits_per_core=10, method=method)
    check(mapping)
    assert len(mapping.assignments) == len(mapping.timeslices)


def test_hqa_needs_fewer_communications_than_the_baselines():
    pairs = random_pairs(60, 60, seed=2)
    counts = {m: map_to_cores(pairs, 6, 10, method=m).communications for m in METHODS}
    assert counts["hqa"] < counts["fgp-roee"] < counts["naive"]
    _, upper = random_circuit_bounds(6, 60, len(pairs))
    assert counts["naive"] <= 1.1 * upper


def test_openqasm_input_and_registers():
    source = """OPENQASM 3.0;
include "stdgates.inc";
qubit[3] a;
qubit[3] b;
h a[0];
cx a[0], b[0];
cx a[1], b[1];
ccx a[2], b[2], a[0];
cx b[0], a[2];
"""
    mapping = map_to_cores(source, cores=2, qubits_per_core=3)
    check(mapping)
    assert [q.label() for q in mapping.qubits] == ["a[0]", "a[1]", "a[2]", "b[0]", "b[1]", "b[2]"]


def test_spare_room_and_errors():
    pairs = random_pairs(10, 10, seed=3)
    check(map_to_cores(pairs, cores=3, qubits_per_core=5, num_qubits=10))
    with pytest.raises(ValueError, match="do not fit"):
        map_to_cores(pairs, cores=2, qubits_per_core=4)
    with pytest.raises(ValueError, match="unknown mapping method"):
        map_to_cores(pairs, 4, 4, method="magic")
    with pytest.raises(ValueError, match="at least two qubits"):
        map_to_cores(pairs, 10, 1)


@pytest.mark.parametrize(
    "topology", ["line", "ring", "grid", "grid:2x3", [(0, 1), (1, 2), (1, 3), (3, 4), (3, 5)]]
)
def test_sparse_topologies(topology):
    pairs = random_pairs(36, 20, seed=4)
    for method in METHODS:
        mapping = map_to_cores(pairs, cores=6, qubits_per_core=6, method=method, topology=topology)
        check(mapping)
        assert mapping.communications >= mapping.moved_qubits


def test_counting_links_pays_off_on_a_line():
    pairs = random_pairs(48, 40, seed=5)
    hqa = map_to_cores(pairs, 8, 6, method="hqa", topology="line").communications
    uniform = map_to_cores(pairs, 8, 6, method="hqa-uniform-cost", topology="line").communications
    assert hqa < uniform


def test_topology_errors():
    pairs = random_pairs(8, 4, seed=6)
    with pytest.raises(ValueError, match="do not connect"):
        map_to_cores(pairs, 4, 2, topology=[(0, 1), (2, 3)])
    with pytest.raises(ValueError, match="not 4"):
        map_to_cores(pairs, 4, 2, topology="grid:3x3")
    with pytest.raises(ValueError, match="unknown topology"):
        map_to_cores(pairs, 4, 2, topology="torus")


def simulate(instructions, n):
    """The state after ``instructions`` (name, qubit indices, params) on ``n`` qubits
    from |0…0⟩. A run of teleportations moves states all at once."""
    state = np.zeros((2,) * n, dtype=complex)
    state[(0,) * n] = 1
    i = 0
    while i < len(instructions):
        name, qubits, params = instructions[i]
        if name == "teleport":
            moves = {}
            while i < len(instructions) and instructions[i][0] == "teleport":
                moves[instructions[i][1][1]] = instructions[i][1][0]
                i += 1
            # Axis `to` takes the state of axis `from`; the slots left empty
            # take the empty destinations' |0⟩.
            vacated = [f for f in moves.values() if f not in moves]
            fresh = [t for t in moves if t not in moves.values()]
            source = list(range(n))
            for to, frm in moves.items():
                source[to] = frm
            for v, f in zip(vacated, fresh):
                source[v] = f
            state = np.transpose(state, source)
            continue
        if name != "barrier":
            k = len(qubits)
            gate = gate_matrix(name, params).reshape((2,) * (2 * k))
            state = np.tensordot(gate, state, axes=(list(range(k, 2 * k)), qubits))
            state = np.moveaxis(state, list(range(k)), qubits)
        i += 1
    return state


QFT_LIKE = """OPENQASM 2.0;
include "qelib1.inc";
qreg q[7];
creg c[7];
u3(0.3, 0.2, 0.1) q[0];
u3(1.3, 0.5, 0.1) q[3];
u3(0.7, 0.4, 0.9) q[5];
"""
QFT_LIKE += "".join(
    f"h q[{a}];\n" + "".join(f"cp({0.5 ** (b - a)}) q[{b}], q[{a}];\n" for b in range(a + 1, 7))
    for a in range(7)
)
QFT_LIKE += "ccx q[0], q[3], q[6];\nmeasure q[2] -> c[2];\n"


@pytest.mark.parametrize("chip", ["line", "ring", "grid:2x2"])
@pytest.mark.parametrize("method", ["hqa", "fgp-roee"])
def test_routed_circuits_run_the_same_computation(chip, method):
    result = map_and_route(QFT_LIKE, cores=2, qubits_per_core=4, chip=chip, method=method)
    k = 4
    couplings = set(result.chip)
    physical = []
    for inst in result.circuit.instructions:
        qs = [q.index for q in inst.qubits]
        if len(qs) == 2:
            a, b = qs
            if inst.name == "teleport":
                assert a // k != b // k
            else:
                assert a // k == b // k and (min(a, b) % k, max(a, b) % k) in couplings
        if inst.name != "measure":
            physical.append((inst.name, qs, list(inst.params)))
    names = [inst.name for inst in result.circuit.instructions]
    assert names.count("swap") == result.swaps
    assert names.count("teleport") == result.mapping.moved_qubits
    assert result.teleportations == result.mapping.communications
    (measure,) = [inst for inst in result.circuit.instructions if inst.name == "measure"]
    assert measure.clbits[0].index == 2

    # The original, decomposed the same way, on 7 qubits.
    from qevedo.compiler.frontend.qasm import parse_qasm
    from qevedo.compiler.passes.decompose import DecomposeToNative

    logical = DecomposeToNative().run(parse_qasm(source=QFT_LIKE))
    index = {q: i for i, q in enumerate(result.mapping.qubits)}
    expected = simulate(
        [
            (inst.name, [index[q] for q in inst.qubits], list(inst.params))
            for inst in logical.instructions
            if inst.name != "measure"
        ],
        7,
    )
    got = simulate(physical, 8)
    # Read the logical qubits where they ended; the empty slot is |0⟩.
    order = list(result.final_layout) + [p for p in range(8) if p not in result.final_layout]
    got = np.transpose(got, order)[(Ellipsis, 0)]
    assert abs(abs(np.vdot(expected.ravel(), got.ravel())) - 1) < 1e-9


def test_all_to_all_chips_need_no_swaps():
    pairs_circuit = "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[8];\n" + "".join(
        f"cx q[{a}], q[{b}];\n" for a, b in random_pairs(8, 10, seed=3)
    )
    result = map_and_route(pairs_circuit, cores=2, qubits_per_core=4, chip="all-to-all")
    assert result.swaps == 0
    mapped = map_to_cores(pairs_circuit, cores=2, qubits_per_core=4)
    assert result.teleportations == mapped.communications
    line = map_and_route(pairs_circuit, cores=2, qubits_per_core=4, chip="line")
    assert line.swaps > 0 and line.teleportations == mapped.communications


def test_chip_errors():
    with pytest.raises(ValueError, match="do not connect"):
        map_and_route(QFT_LIKE, cores=2, qubits_per_core=4, chip=[(0, 1), (2, 3)])
    with pytest.raises(ValueError, match="unknown mapping method"):
        map_and_route(QFT_LIKE, cores=2, qubits_per_core=4, method="magic")
    assert topology_links("grid:2x2", 4) == [(0, 1), (0, 2), (1, 3), (2, 3)]


def test_schedules_wait_for_epr_pairs_on_busy_links():
    # Every gate splits the two cores at first: qubits must cross the link.
    pairs = [(0, 4), (1, 5), (2, 6), (3, 7), (0, 1), (2, 3), (4, 5), (6, 7)] * 3
    result = map_and_route(pairs, cores=2, qubits_per_core=4, chip="all-to-all")
    one = result.schedule(link_capacity=1)
    many = result.schedule(link_capacity=8)
    assert len(one.start) == len(result.circuit.instructions)
    assert sum(one.link_hops.values()) == result.teleportations
    assert result.teleportations > 1
    assert one.duration > many.duration
    assert many.epr_wait < one.epr_wait
    with pytest.raises(ValueError, match="channel"):
        result.schedule(link_capacity=0)


@pytest.mark.parametrize("capacity", [1, 4])
def test_hqa_latency_plans_for_the_timing_and_is_no_slower(capacity):
    timing = Timing(link_capacity=capacity)
    pairs = random_pairs(24, 24, seed=5)
    hqa = map_and_route(pairs, cores=4, qubits_per_core=6, chip="line", timing=timing)
    fast = map_and_route(
        pairs, cores=4, qubits_per_core=6, chip="line", method="hqa-latency", timing=timing
    )
    assert fast.timing == hqa.timing == timing
    assert fast.schedule().duration <= hqa.schedule().duration + 1e-9
    # Fields can be changed per call; other timings can be passed whole.
    assert hqa.schedule(link_capacity=capacity + 1).duration <= hqa.schedule().duration
    assert hqa.schedule(Timing(teleport=0.0, epr_pair=0.0)).epr_wait == 0
    mapped = map_to_cores(pairs, 4, 6, method="hqa-latency", timing=timing)
    check(mapped)


@pytest.mark.parametrize("ports", ["central", {(0, 1): (3, 0), (1, 2): (3, 0)}])
def test_teleportations_through_ports(ports):
    result = map_and_route(QFT_LIKE, cores=3, qubits_per_core=4, chip="line", topology="line", ports=ports)
    k = 4
    physical = []
    for inst in result.circuit.instructions:
        qs = [q.index for q in inst.qubits]
        if inst.name == "teleport":
            a, b = qs
            assert abs(a // k - b // k) == 1  # one link at a time
            if isinstance(ports, dict):
                link = (min(a // k, b // k), max(a // k, b // k))
                slot = dict(zip(link, ports[link]))
                assert (a % k, b % k) == (slot[a // k], slot[b // k])
        if inst.name != "measure":
            physical.append((inst.name, qs, list(inst.params)))
    assert sum(i.name == "teleport" for i in result.circuit.instructions) == result.teleportations
    from qevedo.compiler.frontend.qasm import parse_qasm
    from qevedo.compiler.passes.decompose import DecomposeToNative

    logical = DecomposeToNative().run(parse_qasm(source=QFT_LIKE))
    index = {q: i for i, q in enumerate(result.mapping.qubits)}
    expected = simulate(
        [
            (inst.name, [index[q] for q in inst.qubits], list(inst.params))
            for inst in logical.instructions
            if inst.name != "measure"
        ],
        7,
    )
    got = simulate(physical, 12)
    order = list(result.final_layout) + [p for p in range(12) if p not in result.final_layout]
    got = np.transpose(got, order)[(Ellipsis,) + (0,) * 5]
    assert abs(abs(np.vdot(expected.ravel(), got.ravel())) - 1) < 1e-9
    with pytest.raises(ValueError, match="no port"):
        map_and_route(QFT_LIKE, cores=3, qubits_per_core=4, topology="line", ports={(0, 1): (0, 0)})
