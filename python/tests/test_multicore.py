"""Mapping onto multi-core architectures (qevedo.compiler.multicore)."""

import random

import pytest

from qevedo.compiler.multicore import METHODS, map_to_cores, random_circuit_bounds
from qevedo.compiler.passes import native

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
