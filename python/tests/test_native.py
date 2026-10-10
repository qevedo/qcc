"""The Rust core (qevedo.compiler._core) against the Python reference passes."""

import numpy as np
import pytest

from qevedo.compiler import CompileOptions, compile_source
from qevedo.compiler.ir import Instruction, Qubit
from qevedo.compiler.passes import native
from qevedo.compiler.passes.base import PassManager
from qevedo.compiler.passes.decompose import Lowering

pytestmark = pytest.mark.skipif(not native.available(), reason="the Rust core is not built")

from test_synthesis import DEVICES, device, params_for, random_program  # noqa: E402
from qevedo.compiler.synthesis.gates import GATE_ARITY  # noqa: E402


def test_default_pipeline_is_the_rust_core():
    assert [p.name for p in PassManager.default().passes] == ["native_pipeline"]


@pytest.mark.parametrize("label", sorted(DEVICES))
def test_each_gate_lowers_with_the_same_two_qubit_count(label):
    gates = DEVICES[label]
    reference = Lowering(gates)
    totals = [0, 0]
    for name, n in sorted(GATE_ARITY.items()):
        params = params_for(name)
        python = reference.lower(
            Instruction(name, tuple(Qubit("q", i) for i in range(n)), tuple(params))
        )
        rust = native._core.run_passes(
            [(name, list(range(n)), params, [])], sorted(gates), ["decompose"]
        )
        assert sum(len(q) > 1 for _, q, _, _ in rust) == sum(len(i.qubits) > 1 for i in python), (
            name
        )
        totals[0] += len(rust)
        totals[1] += len(python)
    # Both search the same space, but their single-qubit searches end in
    # different local optima (and the Python one moves with NumPy's LAPACK), so
    # only the sum is compared.
    assert totals[0] <= 1.05 * totals[1]


@pytest.mark.parametrize("label", sorted(DEVICES))
def test_circuits_compile_to_similar_counts(label):
    rng = np.random.default_rng(42)
    totals = np.zeros(4, int)
    for _ in range(10):
        source = random_program(rng, n=4, length=40)
        dev = device(DEVICES[label])
        rust = compile_source(source, CompileOptions(device=dev))
        python = compile_source(
            source, CompileOptions(device=dev, pass_manager=PassManager.reference(dev))
        )
        totals += [
            rust.two_qubit_gate_count(),
            python.two_qubit_gate_count(),
            len(rust.instructions),
            len(python.instructions),
        ]
    assert totals[0] <= 1.02 * totals[1]
    assert totals[2] <= 1.05 * totals[3]


def test_lowering_errors_come_back_as_lowering_errors():
    from qevedo.compiler import LoweringError

    with pytest.raises(LoweringError, match="cannot lower gate 'g'"):
        compile_source('OPENQASM 2.0; include "qelib1.inc"; opaque g a; qreg q[1]; g q[0];')
    with pytest.raises(LoweringError, match="no supported two-qubit gate"):
        compile_source(
            'OPENQASM 3.0; include "stdgates.inc"; qubit[2] q; swap q[0], q[1];',
            CompileOptions(device=device(["rz", "sx", "iswap"])),
        )


def test_registers_and_measurements_round_trip():
    source = """OPENQASM 3.0;
include "stdgates.inc";
qubit[2] a;
qubit[1] b;
bit[3] c;
h a[1];
cx a[1], b[0];
barrier a[0], a[1], b[0];
c[2] = measure b[0];
reset a[0];
"""
    circuit = compile_source(source)
    names = [
        (i.name, [q.label() for q in i.qubits], [x.label() for x in i.clbits])
        for i in circuit.instructions
    ]
    assert ("cx", ["a[1]", "b[0]"], []) in names
    assert ("measure", ["b[0]"], ["c[2]"]) in names
    assert ("reset", ["a[0]"], []) in names
    assert any(n == "barrier" and len(q) == 3 for n, q, _ in names)
