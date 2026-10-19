"""SWAP elision: SWAP gates become a relabelling of the qubits."""

import pytest

from qevedo.compiler import CompileOptions, compile_circuit
from qevedo.compiler.device import DeviceSpec
from qevedo.compiler.frontend.qasm import parse_qasm
from qevedo.compiler.ir import Qubit
from qevedo.compiler.passes import ElideSwaps, PassManager

GATES = {"rz", "sx", "x", "cx"}
PROGRAM = """OPENQASM 2.0;
include "qelib1.inc";
qreg q[3];
creg c[3];
h q[0];
swap q[0], q[1];
cx q[1], q[2];
swap q[1], q[2];
t q[2];
"""
MEASURED = PROGRAM + "measure q -> c;\n"


def device(topology="all_to_all"):
    return DeviceSpec(name="test", topology=topology, qubits=3, native_gates=set(GATES))


def q(i):
    return Qubit("q", i)


def test_swaps_go_when_their_qubits_are_measured_last():
    out = ElideSwaps().run(parse_qasm(source=MEASURED))
    assert [i.name for i in out.instructions] == ["h", "cx", "t", "measure", "measure", "measure"]
    # What the input leaves on q[0] ends on q[1], q[1]'s on q[2], q[2]'s on q[0].
    assert out.output_permutation == {q(0): q(1), q(1): q(2), q(2): q(0)}
    measures = {i.clbits[0].index: i.qubits[0].index for i in out.instructions[3:]}
    assert measures == {0: 1, 1: 2, 2: 0}


def test_swaps_stay_unless_the_outputs_may_be_relabelled():
    circuit = parse_qasm(source=PROGRAM)
    assert ElideSwaps().run(circuit) is circuit
    assert ElideSwaps(relabel_outputs=True).run(circuit).output_permutation
    measured = parse_qasm(source=MEASURED)
    assert ElideSwaps().run(measured, device("linear_chain")) is measured


def test_compiled_circuits_agree_up_to_the_permutation():
    qiskit = pytest.importorskip("qiskit")
    from qiskit.quantum_info import Operator

    from qevedo.compiler import emit_qasm

    circuit = parse_qasm(source=PROGRAM)
    kept = compile_circuit(circuit, CompileOptions(device=device()))
    elided = compile_circuit(circuit, CompileOptions(device=device(), relabel_outputs=True))
    assert kept.output_permutation == {}
    assert (elided.two_qubit_gate_count(), kept.two_qubit_gate_count()) == (1, 5)
    # Moving every state back to its own wire gives the input's unitary.
    restored = qiskit.QuantumCircuit.from_qasm_str(emit_qasm(elided))
    wires = {a.index: b.index for a, b in elided.output_permutation.items()}
    for a in range(3):
        b = wires.get(a, a)
        if a != b:
            restored.swap(a, b)
            other = next(x for x, w in wires.items() if w == a)
            wires[other], wires[a] = b, a
    assert Operator(restored).equiv(Operator(qiskit.QuantumCircuit.from_qasm_str(PROGRAM)))


def test_reference_and_default_pipelines_both_elide():
    circuit = parse_qasm(source=MEASURED)
    for manager in (PassManager.default(device()), PassManager.reference(device())):
        out = manager.run(circuit, device())
        assert out.two_qubit_gate_count() == 1


def test_a_cx_each_way_becomes_one_cx_when_measured():
    program = """OPENQASM 2.0;
include "qelib1.inc";
qreg q[2];
creg c[2];
h q[0];
cx q[0], q[1];
cx q[1], q[0];
measure q -> c;
"""
    circuit = parse_qasm(source=program)
    out = compile_circuit(circuit, CompileOptions(device=device()))
    assert out.two_qubit_gate_count() == 1
    assert out.output_permutation == {q(0): q(1), q(1): q(0)}
    unmeasured = parse_qasm(source=program.replace("measure q -> c;\n", ""))
    assert compile_circuit(unmeasured, CompileOptions(device=device())).two_qubit_gate_count() == 2
