"""QCC IR → OpenQASM 2 or 3, printed by the ``openqasm`` package (3.x)."""

from __future__ import annotations

import math
from fractions import Fraction

import openqasm
from openqasm import ast

from qcc.ir import Circuit, Clbit, Instruction, Qubit

__all__ = ["emit_qasm"]

# IR gate names that are spelled differently in OpenQASM 3's stdgates.inc.
_OPENQASM3_NAMES = {"u": "U"}


def emit_qasm(circuit: Circuit, version: str = "2.0", include_barriers: bool = True) -> str:
    """Serialize ``circuit`` as an OpenQASM ``version`` ("2.0" or "3.0") program.

    Gates keep their IR names and come from ``qelib1.inc`` (OpenQASM 2) or
    ``stdgates.inc`` (OpenQASM 3).
    """
    if version not in ("2.0", "3.0"):
        raise ValueError(f"unsupported OpenQASM version {version!r}")
    openqasm3 = version == "3.0"
    statements: list[ast.Statement] = [ast.Include("stdgates.inc" if openqasm3 else "qelib1.inc")]
    for name, size in circuit.qreg_sizes.items():
        statements.append(ast.QubitDeclaration(ast.Identifier(name), ast.IntegerLiteral(size)))
    for name, size in circuit.creg_sizes.items():
        statements.append(
            ast.ClassicalDeclaration(ast.BitType(ast.IntegerLiteral(size)), ast.Identifier(name))
        )
    for inst in circuit.instructions:
        statement = _statement(inst, openqasm3, include_barriers)
        if statement is not None:
            statements.append(statement)
    # The printer uses OpenQASM 2 syntax (qreg, creg, measure ->) for version 2.0.
    return openqasm.dumps(ast.Program(statements, version=version)).rstrip("\n")


def _statement(inst: Instruction, openqasm3: bool, include_barriers: bool) -> ast.Statement | None:
    name = inst.name.lower()
    if name == "barrier":
        return ast.Barrier([_bit(q) for q in inst.qubits]) if include_barriers else None
    if name == "measure":
        measure = ast.MeasureExpression(_bit(inst.qubits[0]))
        return ast.MeasureStatement(measure, _bit(inst.clbits[0]))
    if name == "reset":
        return ast.Reset(_bit(inst.qubits[0]))
    name = _OPENQASM3_NAMES.get(name, inst.name) if openqasm3 else inst.name
    return ast.GateCall(
        [], ast.Identifier(name), [_angle(p) for p in inst.params], [_bit(q) for q in inst.qubits]
    )


def _bit(bit: Qubit | Clbit) -> ast.IndexedIdentifier:
    return ast.IndexedIdentifier(ast.Identifier(bit.register), [[ast.IntegerLiteral(bit.index)]])


def _angle(value: float) -> ast.Expression:
    """``value`` as a readable expression: a small rational multiple of pi when it is one."""
    ratio = Fraction(value / math.pi).limit_denominator(16)
    if ratio and abs(float(ratio) * math.pi - value) < 1e-12:
        # Build "-3 * pi / 4" rather than "-(3 * pi / 4)".
        expression: ast.Expression = ast.Identifier("pi")
        numerator = abs(ratio.numerator)
        if numerator != 1:
            factor: ast.Expression = ast.IntegerLiteral(numerator)
            if ratio < 0:
                factor = ast.UnaryExpression(ast.UnaryOperator.NEGATE, factor)
            expression = ast.BinaryExpression(ast.BinaryOperator.MULTIPLY, factor, expression)
        elif ratio < 0:
            expression = ast.UnaryExpression(ast.UnaryOperator.NEGATE, expression)
        if ratio.denominator != 1:
            expression = ast.BinaryExpression(
                ast.BinaryOperator.DIVIDE, expression, ast.IntegerLiteral(ratio.denominator)
            )
        return expression
    if value < 0:
        return ast.UnaryExpression(ast.UnaryOperator.NEGATE, ast.FloatLiteral(-value))
    return ast.FloatLiteral(value)
