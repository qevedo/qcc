"""QCC as Qiskit's optimization stage (qevedo.compiler.qiskit_plugin)."""

import numpy as np
import pytest

qiskit = pytest.importorskip("qiskit")

from qiskit import QuantumCircuit, transpile  # noqa: E402
from qiskit.circuit import Parameter  # noqa: E402
from qiskit.providers.fake_provider import GenericBackendV2  # noqa: E402
from qiskit.quantum_info import Operator  # noqa: E402
from qiskit.transpiler import CouplingMap  # noqa: E402
from qiskit.transpiler.preset_passmanagers.plugin import list_stage_plugins  # noqa: E402

from qevedo.compiler.passes import native  # noqa: E402

pytestmark = pytest.mark.skipif(not native.available(), reason="the Rust core is not built")


def random_circuit(n, gates, seed):
    rng = np.random.default_rng(seed)
    qc = QuantumCircuit(n)
    for _ in range(gates):
        kind = rng.integers(4)
        q = [int(x) for x in rng.choice(n, 2, replace=False)]
        if kind == 0:
            qc.cx(*q)
        elif kind == 1:
            qc.rzz(float(rng.uniform(-3, 3)), *q)
        elif kind == 2:
            qc.u(*[float(a) for a in rng.uniform(-3, 3, 3)], q[0])
        else:
            qc.h(q[0])
    return qc


def test_the_plugin_is_registered():
    assert "qcc" in list_stage_plugins("optimization")


@pytest.mark.parametrize("seed", range(5))
def test_routed_circuits_stay_equivalent_and_on_the_coupling_map(seed):
    qc = random_circuit(5, 40, seed)
    coupling = CouplingMap.from_line(5)
    out = transpile(
        qc,
        coupling_map=coupling,
        basis_gates=["rz", "sx", "x", "cx"],
        optimization_level=3,
        optimization_method="qcc",
        seed_transpiler=seed,
    )
    assert Operator.from_circuit(out).equiv(Operator(qc))
    edges = {tuple(e) for e in coupling.get_edges()}
    for inst in out.data:
        if inst.operation.num_qubits == 2:
            pair = tuple(out.find_bit(q).index for q in inst.qubits)
            assert pair in edges or pair[::-1] in edges
    assert set(out.count_ops()) <= {"rz", "sx", "x", "cx"}


def test_directed_targets_and_measurements():
    backend = GenericBackendV2(6, basis_gates=["ecr", "rz", "sx", "x"], coupling_map=[[i, i + 1] for i in range(5)], seed=2)
    qc = random_circuit(6, 30, 7)
    qc.measure_all()
    out = transpile(qc, backend, optimization_level=2, optimization_method="qcc", seed_transpiler=1)
    for inst in out.data:
        if inst.operation.name == "ecr":
            qargs = tuple(out.find_bit(q).index for q in inst.qubits)
            assert backend.target.instruction_supported("ecr", qargs)
    measures = [inst for inst in out.data if inst.operation.name == "measure"]
    assert len(measures) == 6
    assert len({out.find_bit(inst.clbits[0]).index for inst in measures}) == 6


def test_circuits_it_cannot_handle_pass_through():
    theta = Parameter("θ")
    qc = QuantumCircuit(2)
    qc.rz(theta, 0)
    qc.cx(0, 1)
    out = transpile(qc, basis_gates=["rz", "sx", "x", "cx"], optimization_method="qcc")
    assert out.parameters == {theta}
