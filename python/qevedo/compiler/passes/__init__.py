"""Compilation passes."""

from qevedo.compiler.passes.base import Pass, PassManager
from qevedo.compiler.passes.decompose import DecomposeToNative, LoweringError
from qevedo.compiler.passes.optimize import (
    CancelAdjacentInverses,
    CommutativeCancellation,
    MergeSingleQubitGates,
)
from qevedo.compiler.passes.native import NativePipeline
from qevedo.compiler.passes.permutation import ElideSwaps
from qevedo.compiler.passes.phase_folding import PhaseFolding
from qevedo.compiler.passes.resynthesize import ResynthesizeTwoQubitBlocks

__all__ = [
    "CancelAdjacentInverses",
    "CommutativeCancellation",
    "DecomposeToNative",
    "ElideSwaps",
    "LoweringError",
    "MergeSingleQubitGates",
    "NativePipeline",
    "Pass",
    "PassManager",
    "PhaseFolding",
    "ResynthesizeTwoQubitBlocks",
]
