"""QCC command-line interface."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from qevedo.compiler.compiler import CompileOptions, compile_file
from qevedo.compiler.device import DeviceSpec, default_device
from qevedo.compiler.emit.qasm import emit_qasm
from qevedo.compiler.frontend.qasm import QasmFrontendError
from qevedo.compiler.passes.decompose import LoweringError


def _profiles_dir() -> Path:
    return Path(__file__).resolve().parent / "profiles"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="qcc", description="Quantum Compiler Collection")
    sub = parser.add_subparsers(dest="command", required=True)

    compile_cmd = sub.add_parser("compile", help="Compile an OpenQASM 2 or 3 circuit")
    compile_cmd.add_argument("input", type=Path, help="Input .qasm file")
    compile_cmd.add_argument("-o", "--output", type=Path, required=True, help="Output .qasm file")
    compile_cmd.add_argument(
        "--device",
        type=Path,
        default=None,
        help="Device profile YAML (default: built-in all-to-all)",
    )
    compile_cmd.add_argument(
        "--qasm-version",
        choices=["2", "3"],
        default="2",
        help="OpenQASM version of the output (default: 2)",
    )
    compile_cmd.add_argument(
        "-I",
        "--include",
        action="append",
        default=[],
        type=Path,
        metavar="DIR",
        help="Directory searched by include statements (repeatable)",
    )
    compile_cmd.add_argument(
        "--stats",
        action="store_true",
        help="Print gate counts and depth to stderr",
    )

    args = parser.parse_args(argv)
    if args.command == "compile":
        return _run_compile(args)
    return 1


def _run_compile(args: argparse.Namespace) -> int:
    device = DeviceSpec.from_yaml(args.device) if args.device else default_device()
    try:
        options = CompileOptions(device=device, include_paths=args.include)
        compiled = compile_file(args.input, options)
    except (QasmFrontendError, LoweringError) as error:
        print(f"qcc: {error}", file=sys.stderr)
        return 1
    output = emit_qasm(compiled, version=f"{args.qasm_version}.0")
    args.output.write_text(output + "\n", encoding="utf-8")
    if args.stats:
        counts = compiled.gate_counts()
        print(f"device: {device.name}", file=sys.stderr)
        print(f"depth: {compiled.depth()}", file=sys.stderr)
        print(f"2q_gates: {compiled.two_qubit_gate_count()}", file=sys.stderr)
        print(f"gates: {counts}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
