"""The default pipeline, run by the Rust core (``qevedo.compiler._core``).

It does what the Python passes in this package do (lowering, single-qubit
merging, commutative cancellation and two-qubit block resynthesis), the same
way, much faster. The Python passes remain as the reference implementation.
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
    """Lower to the device's native gates and optimize, in Rust."""

    name = "native_pipeline"

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
        native = sorted(g.lower() for g in (device or default_device()).native_gates)
        try:
            lowered = _core.compile(raw, native)
        except ValueError as error:
            raise LoweringError(str(error)) from None
        out = circuit.copy()
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
