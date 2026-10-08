"""Compilation passes."""

from qcc.passes.base import Pass, PassManager
from qcc.passes.decompose import DecomposeToNative
from qcc.passes.optimize import CancelAdjacentInverses

__all__ = [
    "Pass",
    "PassManager",
    "DecomposeToNative",
    "CancelAdjacentInverses",
]
