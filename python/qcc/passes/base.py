"""Pass manager infrastructure."""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Iterable, List, Optional

from qcc.device import DeviceSpec
from qcc.ir import Circuit


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
    def default(cls, device: Optional[DeviceSpec] = None) -> PassManager:
        from qcc.passes.decompose import DecomposeToNative
        from qcc.passes.optimize import CancelAdjacentInverses

        del device  # routing passes will use this later
        return cls([DecomposeToNative(), CancelAdjacentInverses()])
