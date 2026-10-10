"""Circuit-level optimizations: block resynthesis and commutative cancellation."""

import math

import numpy as np
import pytest

from qevedo.compiler.device import DeviceSpec
from qevedo.compiler.ir import Circuit, Instruction, Qubit
from qevedo.compiler.passes import CommutativeCancellation, ResynthesizeTwoQubitBlocks
from qevedo.compiler.passes.decompose import cost
from qevedo.compiler.synthesis.gates import embed, gate_matrix

N = 3
Q = [Qubit("q", i) for i in range(N)]
NATIVE = {
    "ibm": ["rz", "sx", "x", "cx"],
    "cz": ["rz", "rx", "cz"],
    "ions": ["rx", "ry", "rxx"],
}


def device(native):
    return DeviceSpec(name="test", topology="all_to_all", qubits=N, native_gates=set(native))


def circuit(*instructions):
    return Circuit(num_qubits=N, qreg_sizes={"q": N}, instructions=list(instructions))


def gate(name, *qubits, params=()):
    return Instruction(name, tuple(Q[q] for q in qubits), tuple(params))


def unitary(c):
    m = np.eye(2**N, dtype=complex)
    for inst in c.instructions:
        if inst.name == "barrier":
            continue
        m = embed(gate_matrix(inst.name, inst.params), [q.index for q in inst.qubits], N) @ m
    return m


def equal_up_to_phase(a, b):
    k = np.argmax(abs(b))
    phase = a.flat[k] / b.flat[k]
    return abs(abs(phase) - 1) < 1e-8 and np.allclose(a, phase * b, atol=1e-8)


def random_native_circuit(rng, native, length):
    two = next(g for g in native if g in ("cx", "cz", "rxx"))
    singles = [g for g in native if g != two]
    out = []
    for _ in range(length):
        if rng.uniform() < 0.4:
            a, b = (int(x) for x in rng.permutation(N)[:2])
            params = (float(rng.uniform(-3, 3)),) if two == "rxx" else ()
            out.append(gate(two, a, b, params=params))
        else:
            name = singles[rng.integers(len(singles))]
            params = (
                ()
                if name in ("sx", "x")
                else (float(rng.choice([rng.uniform(-3, 3), math.pi / 2])),)
            )
            out.append(gate(name, int(rng.integers(N)), params=params))
    return circuit(*out)


@pytest.mark.parametrize("label", sorted(NATIVE))
@pytest.mark.parametrize("make", [ResynthesizeTwoQubitBlocks, CommutativeCancellation])
def test_passes_preserve_the_unitary_and_never_cost_more(label, make):
    rng = np.random.default_rng(11)
    dev = device(NATIVE[label])
    for _ in range(150):
        c = random_native_circuit(rng, NATIVE[label], int(rng.integers(2, 14)))
        out = make().run(c, dev)
        assert equal_up_to_phase(unitary(out), unitary(c))
        # Fewer two-qubit gates first, then fewer gates.
        assert cost(out.instructions) <= cost(c.instructions)
        assert {i.name for i in out.instructions} <= set(NATIVE[label])


def test_resynthesis_collapses_blocks():
    dev = device(NATIVE["ibm"])
    resynthesize = ResynthesizeTwoQubitBlocks()
    # Two SWAPs on the same pair are the identity: six CX become nothing.
    swaps = [gate("cx", 0, 1), gate("cx", 1, 0), gate("cx", 0, 1)] * 2
    assert resynthesize.run(circuit(*swaps), dev).instructions == []
    # A block with four CXs that is a CZ up to single-qubit gates needs one.
    block = circuit(
        gate("cx", 0, 1), gate("rz", 1, params=(0.3,)), gate("cx", 0, 1),
        gate("cx", 1, 0), gate("rz", 0, params=(0.7,)), gate("cx", 1, 0),
    )  # fmt: skip
    out = resynthesize.run(block, dev)
    assert sum(i.name == "cx" for i in out.instructions) == 2
    assert equal_up_to_phase(unitary(out), unitary(block))
    # Gates on another qubit in between do not split a block.
    split = circuit(gate("cx", 0, 1), gate("sx", 2), gate("cx", 0, 1))
    assert [i.name for i in resynthesize.run(split, dev).instructions] == ["sx"]


def test_commutative_cancellation():
    run = lambda *insts: [  # noqa: E731
        (i.name, [q.index for q in i.qubits], [round(p, 6) for p in i.params])
        for i in CommutativeCancellation().run(circuit(*insts)).instructions
    ]
    # rz commutes with a CX control and merges with the next rz.
    assert run(gate("rz", 0, params=(0.25,)), gate("cx", 0, 1), gate("rz", 0, params=(0.5,))) == [
        ("rz", [0], [0.75]),
        ("cx", [0, 1], []),
    ]
    # ...but not with a CX target.
    assert (
        len(run(gate("rz", 1, params=(0.25,)), gate("cx", 0, 1), gate("rz", 1, params=(0.5,)))) == 3
    )
    # Two CXs cancel through gates that commute with both of them.
    assert run(
        gate("cx", 0, 1),
        gate("rz", 0, params=(0.4,)),
        gate("x", 1),
        gate("cx", 0, 2),
        gate("cx", 0, 1),
    ) == [("rz", [0], [0.4]), ("x", [1], []), ("cx", [0, 2], [])]
    # A gate that commutes on one qubit but not the other blocks the pair.
    assert len(run(gate("cx", 0, 1), gate("cx", 1, 2), gate("cx", 0, 1))) == 3
    # Rotations that add up to nothing disappear; barriers block everything.
    assert run(gate("rz", 0, params=(0.3,)), gate("cz", 0, 1), gate("rz", 0, params=(-0.3,))) == [
        ("cz", [0, 1], [])
    ]
    assert len(run(gate("x", 0), Instruction("barrier", (Q[0],)), gate("x", 0))) == 3
