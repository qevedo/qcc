"""QCC inside Qiskit: a transpiler plugin for the optimization stage.

With ``qevedo-compiler`` and Qiskit installed, Qiskit finds the plugin by
itself::

    from qiskit import transpile

    compiled = transpile(circuit, backend, optimization_method="qcc")

Qiskit still lays out, routes and translates the circuit to the backend's
gates; QCC then optimizes it, with the Rust core's pipeline: phase folding,
single-qubit merging, commutative cancellation and two-qubit block
resynthesis. QCC only
rewrites gates on pairs of qubits that already interact, so the routing
stays valid.

The pass hands the circuit back unchanged when it cannot help safely: control
flow, unbound parameters, gates QCC does not know, or a result that would use
a pair of qubits, or a direction of a two-qubit gate, that the input did not
use and the target does not allow.
"""

from __future__ import annotations

from qiskit.circuit import Barrier, Measure, Reset
from qiskit.circuit.library import get_standard_gate_name_mapping
from qiskit.dagcircuit import DAGCircuit
from qiskit.transpiler import PassManager, TransformationPass
from qiskit.transpiler.preset_passmanagers.plugin import PassManagerStagePlugin

__all__ = ["QCCOptimization", "QCCPlugin"]

# Gates the QCC core reads, by Qiskit's names.
_KNOWN = {
    "id", "x", "y", "z", "h", "s", "sdg", "t", "tdg", "sx", "sxdg",
    "rx", "ry", "rz", "p", "u1", "u2", "u3", "u",
    "cx", "cy", "cz", "ch", "csx", "swap", "iswap", "dcx", "ecr",
    "cp", "cu1", "crx", "cry", "crz", "rxx", "ryy", "rzz", "rzx", "cu3", "cu",
    "ccx", "cswap",
}  # fmt: skip
_NON_GATES = {"measure", "reset", "barrier"}


class QCCOptimization(TransformationPass):
    """Optimize a circuit already in ``basis_gates`` with QCC's Rust pipeline."""

    def __init__(self, basis_gates, target=None):
        super().__init__()
        self.basis_gates = {g for g in basis_gates or () if g in _KNOWN}
        self.target = target

    def run(self, dag: DAGCircuit) -> DAGCircuit:
        from qevedo.compiler import _core

        if not self.basis_gates:
            return dag
        qubit = {q: i for i, q in enumerate(dag.qubits)}
        clbit = {c: i for i, c in enumerate(dag.clbits)}
        raw = []
        for node in dag.topological_op_nodes():
            op = node.op
            if op.name not in _KNOWN | _NON_GATES or getattr(op, "condition", None) is not None:
                return dag
            try:
                params = [float(p) for p in op.params]
            except TypeError:  # an unbound parameter
                return dag
            raw.append(
                (op.name, [qubit[q] for q in node.qargs], params, [clbit[c] for c in node.cargs])
            )
        try:
            out = _core.compile(raw, sorted(self.basis_gates))
        except ValueError:
            return dag
        if not self._allowed(raw, out):
            return dag
        result = dag.copy_empty_like()
        gates = get_standard_gate_name_mapping()
        for name, qs, params, cs in out:
            if name == "measure":
                op = Measure()
            elif name == "reset":
                op = Reset()
            elif name == "barrier":
                op = Barrier(len(qs))
            else:
                op = type(gates[name])(*params)
            result.apply_operation_back(
                op, qargs=tuple(dag.qubits[i] for i in qs), cargs=tuple(dag.clbits[i] for i in cs)
            )
        return result

    def _allowed(self, raw, out) -> bool:
        """Whether every two-qubit gate of ``out`` is on a pair (and in a
        direction) the input used, or the target allows."""
        used = {tuple(qs) for name, qs, _, _ in raw if len(qs) == 2 and name not in _NON_GATES}
        for name, qs, _, _ in out:
            if len(qs) != 2 or name in _NON_GATES or tuple(qs) in used:
                continue
            if self.target is None:
                # Without a target, any direction of a pair already in use.
                if tuple(reversed(qs)) not in used:
                    return False
            elif not self.target.instruction_supported(name, tuple(qs)):
                return False
        return True


class QCCPlugin(PassManagerStagePlugin):
    """The ``"qcc"`` optimization method for Qiskit's ``transpile``."""

    def pass_manager(self, pass_manager_config, optimization_level=None) -> PassManager:
        target = pass_manager_config.target
        basis = pass_manager_config.basis_gates
        if basis is None and target is not None:
            basis = target.operation_names
        return PassManager([QCCOptimization(basis, target)])
