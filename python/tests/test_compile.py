from pathlib import Path

import pytest

from qcc.device import DeviceSpec
from qcc.frontend.qasm import parse_qasm
from qcc.passes.decompose import DecomposeToNative
from qcc.passes.optimize import CancelAdjacentInverses


BELL = Path(__file__).resolve().parents[1] / "examples" / "bell.qasm"
PROFILES = Path(__file__).resolve().parents[1] / "profiles"


def test_parse_bell_circuit():
    circuit = parse_qasm(path=BELL)
    assert circuit.num_qubits == 2
    assert circuit.num_clbits == 2
    assert len(circuit.instructions) == 4
    assert circuit.instructions[0].name == "h"
    assert circuit.instructions[1].name == "cx"


def test_decompose_removes_h():
    circuit = parse_qasm(path=BELL)
    decomposed = DecomposeToNative().run(circuit)
    names = [inst.name for inst in decomposed.instructions]
    assert "h" not in names
    assert "rz" in names
    assert "sx" in names


def test_cancel_double_x():
    from qcc.ir import Circuit, Instruction, Qubit

    circuit = Circuit(num_qubits=1, qreg_sizes={"q": 1})
    q = Qubit("q", 0)
    circuit.extend([Instruction("x", (q,)), Instruction("x", (q,))])
    optimized = CancelAdjacentInverses().run(circuit)
    assert optimized.instructions == []


def test_load_grid_profile():
    spec = DeviceSpec.from_yaml(PROFILES / "grid_2x4.yaml")
    assert spec.qubits == 8
    assert "cx" in spec.native_gates
    assert spec.supports_gate("cx", (0, 1))
    assert not spec.supports_gate("cx", (0, 2))
