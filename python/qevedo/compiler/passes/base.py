"""Pass manager infrastructure."""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Iterable, List, Optional

from qevedo.compiler.device import DeviceSpec
from qevedo.compiler.ir import Circuit


class Pass(ABC):
    """One transformation over a circuit."""

    name: str = "pass"

    @abstractmethod
    def run(self, circuit: Circuit, device: Optional[DeviceSpec] = None) -> Circuit:
        raise NotImplementedError


class PassManager:
    """Run an ordered list of passes."""

    def __init__(self, passes: Iterable[Pass]):
        self.passes: List[Pass] = list(passes)

    def run(self, circuit: Circuit, device: Optional[DeviceSpec] = None) -> Circuit:
        current = circuit.copy()
        for pas in self.passes:
            current = pas.run(current, device)
        return current

    @classmethod
    def default(
        cls, device: Optional[DeviceSpec] = None, relabel_outputs: bool = False
    ) -> PassManager:
        """The default pipeline: the Rust core when it is installed, else the
        Python passes. SWAP gates are elided first; with ``relabel_outputs``,
        even when that leaves the outputs on other wires (see
        ``Circuit.output_permutation``)."""
        from qevedo.compiler.passes.native import NativePipeline, available
        from qevedo.compiler.passes.permutation import ElideSwaps

        if available():
            return cls([ElideSwaps(relabel_outputs), NativePipeline(relabel_outputs)])
        return cls.reference(device, relabel_outputs)

    @classmethod
    def reference(
        cls, device: Optional[DeviceSpec] = None, relabel_outputs: bool = False
    ) -> PassManager:
        """The default pipeline as Python passes, the reference for the Rust core."""
        from qevedo.compiler.passes.decompose import DecomposeToNative
        from qevedo.compiler.passes.optimize import CommutativeCancellation, MergeSingleQubitGates
        from qevedo.compiler.passes.permutation import ElideSwaps
        from qevedo.compiler.passes.phase_folding import PhaseFolding
        from qevedo.compiler.passes.resynthesize import ResynthesizeTwoQubitBlocks

        del device  # routing passes will use this later
        # Each optimization can expose work for the others (merging single-qubit
        # runs brings inverse pairs together, resynthesis joins blocks), so the
        # round runs twice, the second time without resynthesis, as in the Rust
        # core. Phase folding goes first, while rotations are still Z
        # rotations between CNOTs.
        cleanup = [PhaseFolding(), MergeSingleQubitGates(), CommutativeCancellation()]
        return cls(
            [
                ElideSwaps(relabel_outputs),
                DecomposeToNative(),
                *cleanup,
                ResynthesizeTwoQubitBlocks(),
                *cleanup,
                MergeSingleQubitGates(),
            ]
        )
