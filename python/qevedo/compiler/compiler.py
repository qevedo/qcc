"""End-to-end compilation driver."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path

from qevedo.compiler.device import DeviceSpec, default_device
from qevedo.compiler.emit.qasm import emit_qasm
from qevedo.compiler.frontend.qasm import parse_qasm
from qevedo.compiler.ir import Circuit
from qevedo.compiler.passes.base import PassManager


@dataclass
class CompileOptions:
    device: DeviceSpec | None = None
    pass_manager: PassManager | None = None
    #: Directories searched by ``include`` statements, after the input's directory.
    include_paths: Sequence[str | Path] = field(default_factory=tuple)


def compile_circuit(circuit: Circuit, options: CompileOptions | None = None) -> Circuit:
    opts = options or CompileOptions()
    device = opts.device or default_device()
    manager = opts.pass_manager or PassManager.default(device)
    return manager.run(circuit, device)


def compile_file(input_path: str | Path, options: CompileOptions | None = None) -> Circuit:
    """Compile an OpenQASM 2 or 3 file."""
    opts = options or CompileOptions()
    circuit = parse_qasm(path=input_path, include_paths=opts.include_paths)
    return compile_circuit(circuit, opts)


def compile_source(source: str, options: CompileOptions | None = None) -> Circuit:
    """Compile OpenQASM 2 or 3 source text."""
    opts = options or CompileOptions()
    circuit = parse_qasm(source=source, include_paths=opts.include_paths)
    return compile_circuit(circuit, opts)


def compile_to_qasm(
    input_path: str | Path,
    options: CompileOptions | None = None,
    version: str = "2.0",
) -> str:
    """Compile a file and return the result as OpenQASM ``version`` source."""
    return emit_qasm(compile_file(input_path, options), version=version)
