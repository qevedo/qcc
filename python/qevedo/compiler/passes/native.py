"""The default pipeline, run by the Rust core (``qevedo.compiler._core``).

It does what the Python passes in this package do (lowering, phase folding,
single-qubit merging, commutative cancellation and two-qubit block
resynthesis), the same way, much faster. The Python passes remain as the
reference implementation.

Where qubits may be relabelled (see ``ElideSwaps``), the core's block
resynthesis is also permutation-aware: it may synthesize a block U as SWAP·U
and elide the SWAP, when that is cheaper. The Python passes do not do this.
"""

from __future__ import annotations

import dataclasses

from qevedo.compiler.device import DeviceSpec, default_device
from qevedo.compiler.ir import Circuit, Clbit, Instruction, Qubit
from qevedo.compiler.passes.base import Pass
from qevedo.compiler.passes.decompose import LoweringError

try:
    from qevedo.compiler import _core
except ImportError:  # pragma: no cover - only without the compiled extension
    _core = None

__all__ = ["NativePipeline", "available"]


def available() -> bool:
    """Whether the compiled Rust core is installed."""
    return _core is not None


class NativePipeline(Pass):
    """Lower to the device's native gates and optimize, in Rust.

    With ``relabel_outputs``, qubits' states may end on other wires (recorded
    in ``output_permutation``); without it, only when every moved state is
    measured last. For a device that is not all-to-all and lists its
    couplings, the circuit is laid out and routed onto its physical qubits
    (``initial_layout`` and ``output_permutation`` say where each logical
    qubit starts and ends).
    """

    name = "native_pipeline"

    def __init__(self, relabel_outputs: bool = False):
        self.relabel_outputs = relabel_outputs

    def run(self, circuit: Circuit, device: DeviceSpec | None = None) -> Circuit:
        if _core is None:
            raise RuntimeError("the qevedo.compiler._core extension is not installed")
        qubits = _Numbering(circuit.qreg_sizes, Qubit)
        clbits = _Numbering(circuit.creg_sizes, Clbit)
        qubit, clbit = qubits.index, clbits.index
        raw = [
            (
                inst.name,
                [qubit(q) for q in inst.qubits],
                list(inst.params),
                [clbit(c) for c in inst.clbits],
            )
            for inst in circuit.instructions
        ]
        device = device or default_device()
        native = sorted(g.lower() for g in device.native_gates)
        if device.topology != "all_to_all" and device.edges:
            return self._on_device(circuit, raw, native, device, qubits, clbits)
        if device.topology != "all_to_all":
            relabel = "never"
        else:
            relabel = "always" if self.relabel_outputs else "if_measured"
        try:
            lowered, wires = _core.compile_relabelling(raw, native, relabel)
        except ValueError as error:
            raise LoweringError(str(error)) from None
        # Compose with where earlier passes left each state: every qubit
        # the circuit uses has been numbered by now.
        moved = {}
        for q in list(qubits.bits) + list(circuit.output_permutation):
            w = qubit(circuit.output_permutation.get(q, q))
            final = qubits.bits[wires[w] if w < len(wires) else w]
            if final != q:
                moved[q] = final
        q_bits, c_bits = qubits.bits, clbits.bits
        instructions = [
            Instruction(
                name,
                tuple([q_bits[i] for i in qs]),
                tuple(params),
                tuple([c_bits[i] for i in cs]),
            )
            for name, qs, params, cs in lowered
        ]
        return dataclasses.replace(
            circuit,
            qreg_sizes=dict(circuit.qreg_sizes),
            creg_sizes=dict(circuit.creg_sizes),
            instructions=instructions,
            output_permutation=moved,
        )


    @staticmethod
    def _on_device(circuit, raw, native, device, qubits, clbits) -> Circuit:
        """Lay out and route on the device's coupled pairs, on its physical
        qubits: a register ``q`` of ``device.qubits``."""
        try:
            lowered, initial, final = _core.compile_for_device(
                raw, native, device.qubits, list(device.edges), device.directed_edges
            )
        except ValueError as error:
            raise LoweringError(str(error)) from None
        physical = [Qubit("q", i) for i in range(device.qubits)]
        c_bits = clbits.bits
        instructions = [
            Instruction(
                name,
                tuple([physical[i] for i in qs]),
                tuple(params),
                tuple([c_bits[i] for i in cs]),
            )
            for name, qs, params, cs in lowered
        ]
        logical = qubits.bits[: len(initial)]
        return dataclasses.replace(
            circuit,
            num_qubits=device.qubits,
            qreg_sizes={"q": device.qubits},
            creg_sizes=dict(circuit.creg_sizes),
            instructions=instructions,
            initial_layout={q: physical[p] for q, p in zip(logical, initial)},
            output_permutation={q: physical[p] for q, p in zip(logical, final)},
        )


class _Numbering:
    """Numbers the bits of named registers, in declaration order; bits outside
    them follow, in the order they are met."""

    def __init__(self, sizes: dict[str, int], kind: type):
        self.bits: list[Qubit | Clbit] = [
            kind(name, i) for name, size in sizes.items() for i in range(size)
        ]
        self.indices = {bit: i for i, bit in enumerate(self.bits)}

    def index(self, bit: Qubit | Clbit) -> int:
        i = self.indices.get(bit)
        if i is None:
            i = self.indices[bit] = len(self.bits)
            self.bits.append(bit)
        return i
