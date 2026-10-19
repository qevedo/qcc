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
    measured last. Never for devices that are not all-to-all.
    """

    name = "native_pipeline"

    def __init__(self, relabel_outputs: bool = False):
        self.relabel_outputs = relabel_outputs

    def run(self, circuit: Circuit, device: DeviceSpec | None = None) -> Circuit:
        if _core is None:
            raise RuntimeError("the qevedo.compiler._core extension is not installed")
        qubits = _Numbering(circuit.qreg_sizes)
        clbits = _Numbering(circuit.creg_sizes)
        raw = [
            (
                inst.name,
                [qubits.index(q) for q in inst.qubits],
                list(inst.params),
                [clbits.index(c) for c in inst.clbits],
            )
            for inst in circuit.instructions
        ]
        device = device or default_device()
        native = sorted(g.lower() for g in device.native_gates)
        if device.topology != "all_to_all":
            relabel = "never"
        else:
            relabel = "always" if self.relabel_outputs else "if_measured"
        try:
            lowered, wires = _core.compile_relabelling(raw, native, relabel)
        except ValueError as error:
            raise LoweringError(str(error)) from None
        out = circuit.copy()
        # Compose with where earlier passes left each state.
        moved = {}
        for q in {q for inst in circuit.instructions for q in inst.qubits} | set(
            circuit.output_permutation
        ):
            w = qubits.index(circuit.output_permutation.get(q, q))
            final = qubits.bit(wires[w] if w < len(wires) else w, Qubit)
            if final != q:
                moved[q] = final
        out.output_permutation = moved
        out.instructions = [
            Instruction(
                name,
                tuple(qubits.bit(i, Qubit) for i in qs),
                tuple(params),
                tuple(clbits.bit(i, Clbit) for i in cs),
            )
            for name, qs, params, cs in lowered
        ]
        return out


class _Numbering:
    """Numbers the bits of named registers, in declaration order."""

    def __init__(self, sizes: dict[str, int]):
        self.bits: list[tuple[str, int]] = [
            (name, i) for name, size in sizes.items() for i in range(size)
        ]
        self.indices = {bit: i for i, bit in enumerate(self.bits)}

    def index(self, bit: Qubit | Clbit) -> int:
        key = (bit.register, bit.index)
        if key not in self.indices:
            self.indices[key] = len(self.bits)
            self.bits.append(key)
        return self.indices[key]

    def bit(self, index: int, kind: type) -> Qubit | Clbit:
        register, i = self.bits[index]
        return kind(register, i)
