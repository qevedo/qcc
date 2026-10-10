"""Gate synthesis, checked against the qvd simulator (qevedo-simulator)."""

import math

import numpy as np
import pytest

from qevedo.compiler import CompileOptions, QasmFrontendError, compile_source, emit_qasm
from qevedo.compiler.device import DeviceSpec
from qevedo.compiler.ir import Circuit, Instruction, Qubit
from qevedo.compiler.passes import CancelAdjacentInverses, DecomposeToNative, LoweringError
from qevedo.compiler.passes.decompose import DEFINITIONS
from qevedo.compiler.synthesis.gates import _PARAMETRIC, GATE_ARITY, embed, gate_matrix
from qevedo.compiler.synthesis.one_qubit import bases_for, euler_zyz, ops_matrix, synthesize_1q
from qevedo.compiler.synthesis.two_qubit import (
    _steps_matrix,
    canonical,
    kak,
    synthesize_2q,
    two_qubit_basis,
)

simulator = pytest.importorskip("qevedo.simulator")

QUARTER = math.pi / 4

# Standard definitions of the gates the simulator does not know by name.
QASM_DEFINITIONS = """
gate rzx(t) a, b { h b; cx a, b; rz(t) b; cx a, b; h b; }
gate ecr a, b { rzx(pi/4) a, b; x a; rzx(-pi/4) a, b; }
gate iswap a, b { s a; s b; h a; cx a, b; cx b, a; h b; }
gate dcx a, b { cx a, b; cx b, a; }
gate ryy(t) a, b { rx(pi/2) a; rx(pi/2) b; cx a, b; rz(t) b; cx a, b; rx(-pi/2) a; rx(-pi/2) b; }
"""


def haar(n, rng):
    z = (rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))) / math.sqrt(2)
    q, r = np.linalg.qr(z)
    return q * (np.diag(r) / abs(np.diag(r)))


def equal_up_to_phase(a, b, atol=1e-8):
    k = np.argmax(abs(b))
    phase = a.flat[k] / b.flat[k]
    return abs(abs(phase) - 1) < atol and np.allclose(a, phase * b, atol=atol)


def params_for(name):
    count = _PARAMETRIC[name][0] if name in _PARAMETRIC else 0
    return [0.37 + 0.41 * i for i in range(count)]


def simulator_unitary(name, params, n):
    """The simulator's matrix of gate ``name``, in the big-endian convention of ``gates``."""
    args = f"({', '.join(map(repr, params))})" if params else ""
    # Qubit 0 is the least significant bit of the simulator's index: applying the
    # gate to q[n-1], ..., q[0] gives the big-endian matrix.
    operands = ", ".join(f"q[{n - 1 - i}]" for i in range(n))
    columns = []
    for j in range(2**n):
        flips = "".join(f"x q[{b}];\n" for b in range(n) if j >> b & 1)
        source = (
            f'OPENQASM 2.0;\ninclude "qelib1.inc";\n{QASM_DEFINITIONS}qreg q[{n}];\n'
            f"{flips}{name}{args} {operands};\n"
        )
        columns.append(np.array(simulator.statevector(simulator.Circuit.from_qasm(source))))
    return np.array(columns).T


@pytest.mark.parametrize("name", sorted(GATE_ARITY))
def test_gate_matrices_match_the_simulator(name):
    params = params_for(name)
    expected = simulator_unitary(name, params, GATE_ARITY[name])
    assert equal_up_to_phase(gate_matrix(name, params), expected)


def definition_matrix(name, params):
    n = 3 if name in ("ccx", "cswap", "rccx") else 2
    m = np.eye(2**n, dtype=complex)
    for sub, sub_params, operands in DEFINITIONS[name](*params):
        m = embed(gate_matrix(sub, sub_params), operands, n) @ m
    return m, n


@pytest.mark.parametrize("name", sorted(DEFINITIONS))
def test_definitions(name):
    params = params_for(name) if name in _PARAMETRIC else []
    m, n = definition_matrix(name, params)
    # rccx is not in GATE_ARITY, but the simulator knows it from qelib1.inc.
    assert equal_up_to_phase(m, simulator_unitary(name, params, n))


# -- single qubit -----------------------------------------------------------------

BASES_1Q = {
    "zsx": ["rz", "sx", "x"],
    "zsx without x": ["rz", "sx"],
    "zxz": ["rz", "rx"],
    "zyz": ["rz", "ry"],
    "xyx": ["rx", "ry"],
    "u3": ["u3"],
    "p and sx": ["p", "sx"],
}


def test_euler_angles():
    rng = np.random.default_rng(1)
    for _ in range(200):
        u = haar(2, rng)
        theta, phi, lam = euler_zyz(u)
        assert 0 <= theta <= math.pi
        rebuilt = (
            gate_matrix("rz", (phi,)) @ gate_matrix("ry", (theta,)) @ gate_matrix("rz", (lam,))
        )
        assert equal_up_to_phase(rebuilt, u)


@pytest.mark.parametrize("label", sorted(BASES_1Q))
def test_single_qubit_synthesis_is_exact_and_short(label):
    native = BASES_1Q[label]
    bases = bases_for(native)
    rng = np.random.default_rng(2)
    samples = [haar(2, rng) for _ in range(200)]
    samples += [gate_matrix(g, params_for(g)) for g, n in GATE_ARITY.items() if n == 1]
    for u in samples:
        ops = synthesize_1q(u, bases)
        assert all(name in native for name, _ in ops)
        assert equal_up_to_phase(ops_matrix(ops), u)
        assert len(ops) <= (1 if label == "u3" else 5 if "sx" in native else 3)


def test_single_qubit_special_cases():
    zsx = bases_for(["rz", "sx", "x"])
    assert synthesize_1q(np.eye(2), zsx) == []
    assert [name for name, _ in synthesize_1q(gate_matrix("rz", (0.3,)), zsx)] == ["rz"]
    assert [name for name, _ in synthesize_1q(gate_matrix("h"), zsx)] == ["rz", "sx", "rz"]
    assert [name for name, _ in synthesize_1q(gate_matrix("x"), zsx)] == ["x"]
    assert [name for name, _ in synthesize_1q(gate_matrix("y"), zsx)] == ["x", "rz"]
    assert synthesize_1q(gate_matrix("sx"), zsx) == [("sx", ())]


# -- two qubits -------------------------------------------------------------------


def test_kak_reconstructs_and_is_canonical():
    rng = np.random.default_rng(3)
    for _ in range(300):
        u = haar(4, rng)
        k = kak(u)
        a, b, c = k.point
        assert QUARTER + 1e-9 >= a >= b - 1e-9 and b >= abs(c) - 1e-9
        assert np.allclose(k.matrix(), u, atol=1e-9)


@pytest.mark.parametrize(
    "name, point",
    [
        ("cx", (QUARTER, 0, 0)),
        ("cz", (QUARTER, 0, 0)),
        ("ecr", (QUARTER, 0, 0)),
        ("iswap", (QUARTER, QUARTER, 0)),
        ("dcx", (QUARTER, QUARTER, 0)),
        ("swap", (QUARTER, QUARTER, QUARTER)),
    ],
)
def test_weyl_points_of_standard_gates(name, point):
    assert np.allclose(kak(gate_matrix(name)).point, point, atol=1e-9)


def test_points_are_local_invariants():
    rng = np.random.default_rng(4)
    for _ in range(50):
        point = canonical(0.6, 0.3, -0.1)
        dressed = np.kron(haar(2, rng), haar(2, rng)) @ point @ np.kron(haar(2, rng), haar(2, rng))
        assert np.allclose(kak(dressed).point, (0.6, 0.3, -0.1), atol=1e-9)


@pytest.mark.parametrize("native", ["cx", "cz", "cy", "ch", "ecr", "rzz", "rxx", "ryy", "rzx"])
def test_two_qubit_synthesis_is_exact(native):
    basis = two_qubit_basis(native)
    rng = np.random.default_rng(5)
    samples = [haar(4, rng) for _ in range(100)]
    samples += [gate_matrix(g, params_for(g)) for g, n in GATE_ARITY.items() if n == 2]
    for u in samples:
        steps = synthesize_2q(u, basis)
        assert all(name in ("u", native) for name, _, _ in steps)
        assert equal_up_to_phase(_steps_matrix(steps), u)


@pytest.mark.parametrize(
    "gate, cx_count",
    [
        ("cz", 1),
        ("ch", 1),
        ("ecr", 1),
        ("crz", 2),
        ("cp", 2),
        ("rzz", 2),
        ("iswap", 2),
        ("dcx", 2),
        ("swap", 3),
        ("cu", 2),
    ],
)
def test_cx_counts_are_optimal(gate, cx_count):
    steps = synthesize_2q(gate_matrix(gate, params_for(gate)), two_qubit_basis("cx"))
    assert sum(name == "cx" for name, _, _ in steps) == cx_count


def test_ising_counts():
    rzz = two_qubit_basis("rzz")
    count = lambda u: sum(name == "rzz" for name, _, _ in synthesize_2q(u, rzz))
    assert count(gate_matrix("rzz", (0.4,))) == 1
    assert count(gate_matrix("cx")) == 1
    assert count(gate_matrix("iswap")) == 2
    assert count(gate_matrix("swap")) == 3


def test_unsupported_two_qubit_natives_are_reported():
    assert two_qubit_basis("iswap") is None
    with pytest.raises(LoweringError, match="no supported two-qubit gate"):
        compile_source(
            'OPENQASM 3.0; include "stdgates.inc"; qubit[2] q; swap q[0], q[1];',
            CompileOptions(device=device(["rz", "sx", "iswap"])),
        )


# -- whole circuits ---------------------------------------------------------------


def device(native):
    return DeviceSpec(name="test", topology="all_to_all", qubits=4, native_gates=set(native))


DEVICES = {
    "ibm (rz sx x cx)": ["rz", "sx", "x", "cx"],
    "cz with rx": ["rz", "rx", "cz"],
    "ecr with u3": ["u3", "ecr"],
    "trapped ion (rx ry rxx)": ["rx", "ry", "rxx"],
    "rzz with p and sx": ["p", "sx", "rzz"],
}


def _accepted(name):
    """Whether OpenQASM 2 programs can use ``name`` (qelib1.inc) and the simulator knows it."""
    if name in ("ecr", "iswap", "dcx", "ryy", "rzx"):
        return False
    args = "(" + ", ".join("0.1" for _ in params_for(name)) + ")" if params_for(name) else ""
    operands = ", ".join(f"q[{i}]" for i in range(GATE_ARITY[name]))
    try:
        compile_source(f'OPENQASM 2.0; include "qelib1.inc"; qreg q[3]; {name}{args} {operands};')
    except (QasmFrontendError, LoweringError):
        return False
    return True


INPUT_GATES = [g for g in sorted(GATE_ARITY) if _accepted(g)]


def random_program(rng, n=4, length=40):
    lines = ["OPENQASM 2.0;", 'include "qelib1.inc";', f"qreg q[{n}];"]
    for q in range(n):
        lines.append(f"u3({rng.uniform(0, 3)}, {rng.uniform(-3, 3)}, {rng.uniform(-3, 3)}) q[{q}];")
    for _ in range(length):
        name = INPUT_GATES[rng.integers(len(INPUT_GATES))]
        qubits = rng.permutation(n)[: GATE_ARITY[name]]
        params = ", ".join(str(rng.uniform(-3, 3)) for _ in params_for(name))
        args = f"({params})" if params else ""
        lines.append(f"{name}{args} " + ", ".join(f"q[{q}]" for q in qubits) + ";")
    return "\n".join(lines) + "\n"


def state(source):
    return np.array(simulator.statevector(simulator.Circuit.from_qasm(source)))


@pytest.mark.parametrize("label", sorted(DEVICES))
def test_compiled_circuits_are_equivalent_and_native(label):
    native = DEVICES[label]
    options = CompileOptions(device=device(native))
    rng = np.random.default_rng(6)
    for _ in range(8):
        source = random_program(rng)
        compiled = compile_source(source, options)
        assert {inst.name for inst in compiled.instructions} <= set(native)
        text = emit_qasm(compiled).replace(
            'include "qelib1.inc";', 'include "qelib1.inc";\n' + QASM_DEFINITIONS
        )
        overlap = abs(np.vdot(state(source), state(text)))
        assert overlap == pytest.approx(1.0, abs=1e-8)


def test_gate_counts_after_compilation():
    options = CompileOptions(device=device(["rz", "sx", "x", "cx"]))
    swap = compile_source(
        'OPENQASM 3.0; include "stdgates.inc"; qubit[2] q; swap q[0], q[1];', options
    )
    assert swap.gate_counts() == {"cx": 3}
    bell = compile_source(
        'OPENQASM 3.0; include "stdgates.inc"; qubit[2] q; h q[0]; cx q[0], q[1];', options
    )
    assert bell.gate_counts() == {"rz": 2, "sx": 1, "cx": 1}
    toffoli = compile_source(
        'OPENQASM 3.0; include "stdgates.inc"; qubit[3] q; ccx q[0], q[1], q[2];', options
    )
    assert toffoli.gate_counts()["cx"] == 6


# -- regressions ------------------------------------------------------------------


def test_u_and_ch_are_lowered_correctly():
    # The old pass dropped the +π shifts of U and lowered ch to an h on the control.
    options = CompileOptions(device=device(["rz", "sx", "x", "cx"]))
    for program in (
        "qubit q; U(0.7, 1.1, -0.4) q;",
        "qubit[2] q; h q[0]; ch q[0], q[1];",
    ):
        source = f'OPENQASM 3.0; include "stdgates.inc"; {program}'
        compiled = emit_qasm(compile_source(source, options))
        assert abs(np.vdot(state(source), state(compiled))) == pytest.approx(1.0, abs=1e-9)


def test_opaque_gates_are_rejected():
    source = 'OPENQASM 2.0; include "qelib1.inc"; opaque g a; qreg q[1]; g q[0];'
    with pytest.raises(LoweringError, match="cannot lower gate 'g'"):
        compile_source(source)


def test_defined_and_library_gates_compile_correctly():
    options = CompileOptions(device=device(["rz", "sx", "x", "cx"]))
    source = """OPENQASM 2.0;
include "qelib1.inc";
gate mine(t) a, b { cu1(t) a, b; ry(t / 3) b; }
qreg q[5];
h q[0]; h q[1]; h q[2]; h q[4];
mine(0.7) q[0], q[3];
c3x q[0], q[1], q[2], q[3];
c3sqrtx q[1], q[2], q[3], q[4];
rc3x q[4], q[0], q[1], q[2];
"""
    compiled = emit_qasm(compile_source(source, options))
    assert abs(np.vdot(state(source), state(compiled))) == pytest.approx(1.0, abs=1e-9)


def test_native_gates_are_kept():
    q = Qubit("q", 0)
    circuit = Circuit(num_qubits=1, qreg_sizes={"q": 1}, instructions=[Instruction("h", (q,))])
    kept = DecomposeToNative().run(circuit, device(["h", "rz", "sx", "cx"]))
    assert [inst.name for inst in kept.instructions] == ["h"]


def _circuit(*instructions):
    return Circuit(num_qubits=3, qreg_sizes={"q": 3}, instructions=list(instructions))


def test_cancellation_rules():
    q0, q1, q2 = (Qubit("q", i) for i in range(3))
    run = lambda *insts: [
        i.name for i in CancelAdjacentInverses().run(_circuit(*insts)).instructions
    ]
    # sx·sx is an x, not the identity.
    assert run(Instruction("sx", (q0,)), Instruction("sx", (q0,))) == ["sx", "sx"]
    assert run(Instruction("sx", (q0,)), Instruction("sxdg", (q0,))) == []
    # Gates on other qubits in between do not block a cancellation...
    assert run(
        Instruction("cx", (q0, q1)), Instruction("h", (q2,)), Instruction("cx", (q0, q1))
    ) == ["h"]
    # ...but gates on the same qubits do.
    assert run(
        Instruction("cx", (q0, q1)), Instruction("h", (q1,)), Instruction("cx", (q0, q1))
    ) == ["cx", "h", "cx"]
    # Symmetric gates cancel in either qubit order; cx does not.
    assert run(Instruction("cz", (q0, q1)), Instruction("cz", (q1, q0))) == []
    assert run(Instruction("cx", (q0, q1)), Instruction("cx", (q1, q0))) == ["cx", "cx"]
    # A controlled rotation by 2π is not the identity.
    assert run(
        Instruction("crx", (q0, q1), (math.pi,)), Instruction("crx", (q0, q1), (math.pi,))
    ) == ["crx", "crx"]
    assert run(Instruction("rz", (q0,), (math.pi,)), Instruction("rz", (q0,), (math.pi,))) == []
