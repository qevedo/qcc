"""Compiling for a device whose qubits are not all coupled: layout and routing."""

from pathlib import Path

import numpy as np
import pytest

from qevedo.compiler import CompileOptions, compile_circuit, emit_qasm
from qevedo.compiler.device import DeviceSpec
from qevedo.compiler.frontend.qasm import parse_qasm

PROFILES = Path(__file__).resolve().parents[1] / "qevedo" / "compiler" / "profiles"
PROGRAM = """OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
h q[0];
cx q[0], q[3];
cx q[1], q[2];
cz q[3], q[1];
t q[2];
cx q[2], q[0];
rzz(0.3) q[0], q[1];
"""


def line(n, directed=False):
    return DeviceSpec(
        name="line",
        topology="linear_chain",
        qubits=n,
        edges=[(i, i + 1) for i in range(n - 1)],
        native_gates={"rz", "sx", "x", "cx"},
        directed_edges=directed,
    )


@pytest.mark.parametrize("directed", [False, True])
def test_routed_circuits_reach_the_same_state(directed):
    qiskit = pytest.importorskip("qiskit")
    from qiskit.quantum_info import Statevector

    device = line(5, directed)
    out = compile_circuit(parse_qasm(source=PROGRAM), CompileOptions(device=device))
    allowed = set(device.edges) | (set() if directed else {(b, a) for a, b in device.edges})
    for inst in out.instructions:
        if len(inst.qubits) == 2:
            assert (inst.qubits[0].index, inst.qubits[1].index) in allowed
    assert out.qreg_sizes == {"q": 5}
    assert sorted(q.index for q in out.initial_layout) == [0, 1, 2, 3]
    # From |0…0⟩ the output reaches the input's state, each logical qubit on
    # the physical qubit it ends on.
    psi = Statevector(qiskit.QuantumCircuit.from_qasm_str(PROGRAM)).data
    final = {q.index: p.index for q, p in out.output_permutation.items()}
    expected = np.zeros(2**5, dtype=complex)
    for index, amplitude in enumerate(psi):
        expected[sum(1 << final[v] for v in range(4) if index >> v & 1)] = amplitude
    got = Statevector(qiskit.QuantumCircuit.from_qasm_str(emit_qasm(out))).data
    assert abs(np.vdot(expected, got)) > 1 - 1e-8


def test_a_circuit_that_fits_needs_no_swap():
    # A chain of CXs fits a line: the perfect layout needs no SWAP.
    program = "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[4];\nh q[0];\n" + "".join(
        f"cx q[{a}], q[{b}];\n" for a, b in [(0, 2), (2, 1), (1, 3)]
    )
    out = compile_circuit(parse_qasm(source=program), CompileOptions(device=line(4)))
    assert out.two_qubit_gate_count() == 3


def test_profiles_with_couplings_are_routed():
    device = DeviceSpec.from_yaml(PROFILES / "grid_2x4.yaml")
    out = compile_circuit(parse_qasm(source=PROGRAM), CompileOptions(device=device))
    pairs = {tuple(sorted(e)) for e in device.edges}
    assert all(
        tuple(sorted(q.index for q in inst.qubits)) in pairs
        for inst in out.instructions
        if len(inst.qubits) == 2
    )
