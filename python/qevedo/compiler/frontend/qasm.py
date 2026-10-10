"""OpenQASM 2/3 → QCC IR frontend, built on the ``openqasm`` package (3.x).

Programs are parsed and semantically checked by ``openqasm`` first, so the
lowering below only deals with valid programs. It supports the subset of the
language that the IR can represent: qubit and bit registers, gate calls with
compile-time parameters (with broadcasting over registers), measurements,
resets, barriers, ``let`` aliases of qubits and ``const`` values. Everything
else (control flow, modifiers, classical computation, timing, calibrations)
raises :class:`QasmFrontendError` rather than being silently dropped.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Sequence
from importlib import resources
from os import PathLike
from pathlib import Path

import openqasm
from openqasm import ast
from openqasm.semantic import Analysis, Scope, SymbolKind

from qevedo.compiler.ir import Circuit, Clbit, Instruction, Qubit
from qevedo.compiler.passes.decompose import DEFINITIONS
from qevedo.compiler.synthesis.gates import is_known_gate

__all__ = ["QasmFrontendError", "parse_qasm"]

# Built-in gates of OpenQASM whose IR names are the lower-case qelib1 names.
_GATE_NAMES = {"U": "u", "CX": "cx"}

_FUNCTIONS: dict[str, Callable[..., float]] = {
    "sin": math.sin,
    "cos": math.cos,
    "tan": math.tan,
    "arcsin": math.asin,
    "arccos": math.acos,
    "arctan": math.atan,
    "asin": math.asin,
    "acos": math.acos,
    "atan": math.atan,
    "exp": math.exp,
    "log": math.log,
    "ln": math.log,
    "sqrt": math.sqrt,
    "floor": math.floor,
    "ceiling": math.ceil,
    "pow": math.pow,
    "mod": math.fmod,
}

_BINARY: dict[ast.BinaryOperator, Callable[[float, float], float]] = {
    ast.BinaryOperator.ADD: lambda a, b: a + b,
    ast.BinaryOperator.SUBTRACT: lambda a, b: a - b,
    ast.BinaryOperator.MULTIPLY: lambda a, b: a * b,
    ast.BinaryOperator.DIVIDE: lambda a, b: a / b,
    ast.BinaryOperator.POWER: lambda a, b: a**b,
    ast.BinaryOperator.MODULO: math.fmod,
}


class QasmFrontendError(Exception):
    """The program is invalid, or uses a feature the QCC IR cannot represent."""

    def __init__(self, message: str, node: ast.Node | None = None):
        location = f"{node.span}: " if node is not None and node.span is not None else ""
        super().__init__(location + message)


def parse_qasm(
    source: str | None = None,
    path: str | PathLike[str] | None = None,
    include_paths: Sequence[str | PathLike[str]] = (),
) -> Circuit:
    """Parse OpenQASM 2 or 3 source text, or a file, into a QCC circuit."""
    if (source is None) == (path is None):
        raise ValueError("pass exactly one of source and path")
    filename = None
    if path is not None:
        filename = str(path)
        source = Path(path).read_text(encoding="utf-8")
    assert source is not None
    try:
        program = openqasm.parse(source, filename=filename)
    except openqasm.QasmError as error:
        raise QasmFrontendError(str(error)) from error
    analysis = openqasm.analyze(program, filename=filename, include_paths=include_paths)
    if analysis.errors:
        raise QasmFrontendError("\n".join(str(error) for error in analysis.errors))
    lowering = _Lowering(analysis)
    for include in analysis.includes:
        if include in _LIBRARIES:
            library = resources.files("openqasm.stdlib").joinpath(include).read_text("utf-8")
            lowering.collect(openqasm.parse(library, ignore_version=True), library=True)
        else:
            source_text = Path(include).read_text(encoding="utf-8")
            lowering.collect(openqasm.parse(source_text, filename=include, ignore_version=True))
    lowering.collect(program)
    return lowering.lower(program)


_LIBRARIES = ("qelib1.inc", "stdgates.inc")


class _Lowering:
    def __init__(self, analysis: Analysis):
        self.analysis = analysis
        self.circuit = Circuit()
        self.aliases: dict[str, list[Qubit]] = {}
        # Gates defined by the program and its own includes, which are always
        # expanded, and by the standard libraries, which are expanded only when
        # the compiler has no better way to lower them.
        self.definitions: dict[str, ast.GateDefinition] = {}
        self.library_definitions: dict[str, ast.GateDefinition] = {}
        # Values of the parameters of the gate definitions being expanded.
        self.bindings: list[dict[str, float]] = []

    def collect(self, program: ast.Program, library: bool = False) -> None:
        target = self.library_definitions if library else self.definitions
        for statement in program.statements:
            if isinstance(statement, ast.GateDefinition):
                target[statement.name.name] = statement

    def lower(self, program: ast.Program) -> Circuit:
        for statement in program.statements:
            self.statement(statement)
        return self.circuit

    # -- statements -------------------------------------------------------------

    def statement(self, node: ast.Statement) -> None:
        if isinstance(
            node,
            (
                ast.Include,
                ast.Pragma,
                ast.CalibrationGrammar,
                ast.GateDefinition,
                ast.OpaqueDefinition,
                ast.ConstDeclaration,
                ast.IODeclaration,
                ast.Nop,
            ),
        ):
            # Definitions are known to the analysis; gates stay named
            # instructions for the passes to decompose. Inputs are only an
            # error where a compile-time value is needed.
            return
        if isinstance(node, ast.Block):
            for statement in node.statements:
                self.statement(statement)
        elif isinstance(node, ast.QubitDeclaration):
            self.declare(node.name.name, node.size, self.circuit.qreg_sizes, node)
            self.circuit.num_qubits += self.circuit.qreg_sizes[node.name.name]
        elif isinstance(node, ast.ClassicalDeclaration):
            self.classical_declaration(node)
        elif isinstance(node, ast.AliasDeclaration):
            self.aliases[node.name.name] = self.alias_qubits(node.value)
        elif isinstance(node, ast.GateCall):
            self.gate_call(node)
        elif isinstance(node, ast.GlobalPhase):
            if node.modifiers or node.qubits:
                raise QasmFrontendError("controlled global phases are not supported", node)
            # An uncontrolled global phase has no observable effect.
        elif isinstance(node, ast.MeasureStatement):
            self.measure(node.measure, node.target, node)
        elif isinstance(node, ast.Reset):
            for qubit in self.qubits(node.operand):
                self.circuit.append(Instruction("reset", (qubit,)))
        elif isinstance(node, ast.Barrier):
            self.barrier(node)
        else:
            raise QasmFrontendError(f"{_describe(node)} is not supported by qcc yet", node)

    def declare(
        self, name: str, size: ast.Expression | None, sizes: dict[str, int], node: ast.Node
    ) -> None:
        sizes[name] = 1 if size is None else self.integer(size, node)

    def classical_declaration(self, node: ast.ClassicalDeclaration) -> None:
        if not isinstance(node.type, ast.BitType):
            raise QasmFrontendError("classical variables other than bits are not supported", node)
        name = node.name.name
        self.declare(name, node.type.size, self.circuit.creg_sizes, node)
        self.circuit.num_clbits += self.circuit.creg_sizes[name]
        if isinstance(node.init, ast.MeasureExpression):
            self.measure(node.init, node.name, node)
        elif node.init is not None:
            raise QasmFrontendError("initialised bit registers are not supported", node)

    def gate_call(self, node: ast.GateCall) -> None:
        if node.modifiers:
            raise QasmFrontendError(
                "gate modifiers (inv, pow, ctrl, negctrl) are not supported", node
            )
        if node.duration is not None:
            raise QasmFrontendError("gate durations are not supported", node)
        name = _GATE_NAMES.get(node.name.name, node.name.name)
        params = tuple(self.number(argument) for argument in node.arguments)
        groups = [self.qubits(operand) for operand in node.qubits]
        width = max(len(group) for group in groups)
        for offset in range(width):
            qubits = tuple(group[0] if len(group) == 1 else group[offset] for group in groups)
            self.apply(name, params, qubits, node)

    def apply(
        self, name: str, params: tuple[float, ...], qubits: tuple[Qubit, ...], node: ast.Node
    ) -> None:
        """Append gate ``name``, expanding it when its definition must be used."""
        definition = self.definitions.get(name)
        if definition is None and not _compiler_knows(name):
            definition = self.library_definitions.get(name)
        if definition is None:
            self.circuit.append(Instruction(name, qubits, params))
            return
        bindings = {p.name: value for p, value in zip(definition.parameters, params)}
        targets = {q.name: qubit for q, qubit in zip(definition.qubits, qubits)}
        self.bindings.append(bindings)
        try:
            for statement in definition.body:
                self.body_statement(statement, targets)
        finally:
            self.bindings.pop()

    def body_statement(self, node: ast.Statement, targets: dict[str, Qubit]) -> None:
        if isinstance(node, ast.GateCall):
            if node.modifiers:
                raise QasmFrontendError(
                    "gate modifiers (inv, pow, ctrl, negctrl) are not supported", node
                )
            name = _GATE_NAMES.get(node.name.name, node.name.name)
            params = tuple(self.number(argument) for argument in node.arguments)
            qubits = tuple(targets[operand.name] for operand in node.qubits)  # type: ignore[union-attr]
            self.apply(name, params, qubits, node)
        elif isinstance(node, ast.GlobalPhase):
            if node.modifiers or node.qubits:
                raise QasmFrontendError("controlled global phases are not supported", node)
        elif isinstance(node, ast.Barrier):
            qubits = tuple(targets[operand.name] for operand in node.operands)  # type: ignore[union-attr]
            self.circuit.append(Instruction("barrier", qubits))
        elif not isinstance(node, ast.Nop):
            raise QasmFrontendError(f"{_describe(node)} is not supported in a gate body", node)

    def measure(
        self,
        measure: ast.MeasureExpression | ast.QuantumCall,
        target: ast.Identifier | ast.IndexedIdentifier | None,
        node: ast.Node,
    ) -> None:
        if isinstance(measure, ast.QuantumCall):
            raise QasmFrontendError("measurements through defcal calls are not supported", node)
        if target is None:
            raise QasmFrontendError("measurements must store their result in a bit", node)
        qubits = self.qubits(measure.operand)
        clbits = self.clbits(target)
        if len(qubits) != len(clbits):
            raise QasmFrontendError(f"measuring {len(qubits)} qubits into {len(clbits)} bits", node)
        self.circuit.extend(
            [Instruction("measure", (q,), (), (c,)) for q, c in zip(qubits, clbits, strict=True)]
        )

    def barrier(self, node: ast.Barrier) -> None:
        if not node.operands:
            qubits = [
                Qubit(name, index)
                for name, size in self.circuit.qreg_sizes.items()
                for index in range(size)
            ]
        else:
            qubits = []
            for operand in node.operands:
                for qubit in self.qubits(operand):
                    if qubit not in qubits:
                        qubits.append(qubit)
        self.circuit.append(Instruction("barrier", tuple(qubits)))

    # -- operands -----------------------------------------------------------------

    def qubits(self, operand: ast.Node) -> list[Qubit]:
        if isinstance(operand, ast.HardwareQubit):
            raise QasmFrontendError("hardware qubits are not supported", operand)
        if isinstance(operand, ast.Identifier):
            if operand.name in self.aliases:
                return list(self.aliases[operand.name])
            size = self.circuit.qreg_sizes[operand.name]
            return [Qubit(operand.name, index) for index in range(size)]
        if isinstance(operand, ast.IndexedIdentifier):
            qubits = self.qubits(operand.name)
            for index in operand.indices:
                qubits = self.select(qubits, index, operand)
            return qubits
        raise QasmFrontendError(f"unsupported qubit operand {_describe(operand)}", operand)

    def alias_qubits(self, value: ast.Expression) -> list[Qubit]:
        if isinstance(value, ast.Concatenation):
            return self.alias_qubits(value.lhs) + self.alias_qubits(value.rhs)
        if isinstance(value, ast.IndexExpression):
            return self.select(self.alias_qubits(value.collection), value.index, value)
        if isinstance(value, ast.Identifier) and (
            value.name in self.aliases or value.name in self.circuit.qreg_sizes
        ):
            return self.qubits(value)
        raise QasmFrontendError("only qubit aliases are supported", value)

    def clbits(self, target: ast.Identifier | ast.IndexedIdentifier) -> list[Clbit]:
        name = target.name if isinstance(target, ast.Identifier) else target.name.name
        bits = [Clbit(name, index) for index in range(self.circuit.creg_sizes[name])]
        if isinstance(target, ast.IndexedIdentifier):
            for index in target.indices:
                bits = self.select(bits, index, target)
        return bits

    def select(self, items: list, index: ast.Index, node: ast.Node) -> list:
        """Apply one ``[...]`` to a register given as the list of its elements."""
        if isinstance(index, ast.DiscreteSet):
            return [items[self.position(value, len(items), node)] for value in index.values]
        if len(index) != 1:
            raise QasmFrontendError("multi-dimensional indices are not supported", node)
        element = index[0]
        if not isinstance(element, ast.Range):
            return [items[self.position(element, len(items), node)]]
        size = len(items)
        start = 0 if element.start is None else self.position(element.start, size, node)
        end = size - 1 if element.end is None else self.position(element.end, size, node)
        step = 1 if element.step is None else self.integer(element.step, node)
        if step == 0:
            raise QasmFrontendError("a range step cannot be zero", node)
        stop = end + 1 if step > 0 else end - 1
        return [items[i] for i in range(start, stop, step)]

    def position(self, expression: ast.Expression, size: int, node: ast.Node) -> int:
        value = self.integer(expression, node)
        if value < 0:
            value += size
        if not 0 <= value < size:
            raise QasmFrontendError(f"index {value} is out of range for size {size}", node)
        return value

    # -- compile-time values ----------------------------------------------------------

    def integer(self, expression: ast.Expression, node: ast.Node) -> int:
        value = self.number(expression)
        if value != int(value):
            raise QasmFrontendError("expected an integer", node)
        return int(value)

    def number(self, node: ast.Expression) -> float:
        """Evaluate a gate parameter, size or index known at compile time."""
        if isinstance(node, (ast.IntegerLiteral, ast.FloatLiteral)):
            return node.value
        if isinstance(node, ast.Identifier):
            if self.bindings and node.name in self.bindings[-1]:
                return self.bindings[-1][node.name]
            return self.constant(node)
        if isinstance(node, ast.UnaryExpression) and node.op is ast.UnaryOperator.NEGATE:
            return -self.number(node.operand)
        if isinstance(node, ast.BinaryExpression) and node.op in _BINARY:
            try:
                return _BINARY[node.op](self.number(node.lhs), self.number(node.rhs))
            except (ZeroDivisionError, OverflowError, ValueError) as error:
                raise QasmFrontendError(f"cannot evaluate expression: {error}", node) from error
        if isinstance(node, ast.FunctionCall) and node.name.name in _FUNCTIONS:
            arguments = [self.number(argument) for argument in node.arguments]
            try:
                return float(_FUNCTIONS[node.name.name](*arguments))
            except (ValueError, OverflowError, TypeError) as error:
                raise QasmFrontendError(f"cannot evaluate expression: {error}", node) from error
        if isinstance(node, ast.Cast) and isinstance(
            node.type, (ast.FloatType, ast.AngleType, ast.IntType, ast.UintType)
        ):
            value = self.number(node.operand)
            return (
                float(int(value)) if isinstance(node.type, (ast.IntType, ast.UintType)) else value
            )
        raise QasmFrontendError("value is not a compile-time constant", node)

    def constant(self, identifier: ast.Identifier) -> float:
        scope: Scope | None = self.analysis.globals
        while scope is not None:
            symbol = scope.symbols.get(identifier.name)
            if symbol is not None:
                constant = symbol.kind in (SymbolKind.CONSTANT, SymbolKind.BUILTIN_CONSTANT)
                if constant and isinstance(symbol.value, (int, float)):
                    return float(symbol.value)
                break
            scope = scope.parent
        raise QasmFrontendError(f"'{identifier.name}' is not a compile-time constant", identifier)


def _compiler_knows(name: str) -> bool:
    """Whether the compiler can lower gate ``name`` without its OpenQASM definition."""
    return is_known_gate(name) or name in DEFINITIONS


_DESCRIPTIONS = {
    "If": "'if' statement",
    "For": "'for' loop",
    "While": "'while' loop",
    "Switch": "'switch' statement",
    "Box": "'box' block",
    "Delay": "'delay' instruction",
    "Assignment": "classical assignment",
    "ExpressionStatement": "classical expression",
    "SubroutineDefinition": "subroutine definition",
    "ExternDeclaration": "extern declaration",
    "AliasDeclaration": "alias",
    "CalibrationBlock": "'cal' block",
    "CalibrationDefinition": "'defcal' definition",
    "Return": "'return' statement",
    "End": "'end' statement",
    "Break": "'break' statement",
    "Continue": "'continue' statement",
}


def _describe(node: ast.Node) -> str:
    name = type(node).__name__
    return _DESCRIPTIONS.get(name, name)
