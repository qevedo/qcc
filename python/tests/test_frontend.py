import math

import openqasm
import pytest

from qevedo.compiler.cli import main
from qevedo.compiler.compiler import compile_source
from qevedo.compiler.emit.qasm import emit_qasm
from qevedo.compiler.frontend.qasm import QasmFrontendError, parse_qasm
from qevedo.compiler.ir import Clbit, Instruction, Qubit


def names(circuit):
    return [inst.name for inst in circuit.instructions]


def test_openqasm3_program():
    circuit = parse_qasm(
        source="""
        OPENQASM 3.0;
        include "stdgates.inc";
        const int n = 3;
        qubit[n] q;
        bit[n] c;
        qubit anc;
        h q;
        rz(pi / 4) q[0];
        U(pi, 0, pi) anc;
        reset anc;
        c = measure q;
        """
    )
    assert circuit.qreg_sizes == {"q": 3, "anc": 1}
    assert circuit.num_qubits == 4 and circuit.num_clbits == 3
    assert names(circuit) == ["h", "h", "h", "rz", "u", "reset", "measure", "measure", "measure"]
    assert circuit.instructions[3].params == (math.pi / 4,)
    assert circuit.instructions[6].clbits == (Clbit("c", 0),)


def test_openqasm2_program():
    circuit = parse_qasm(
        source='OPENQASM 2.0;\ninclude "qelib1.inc";\nqreg q[2];\ncreg c[2];\n'
        "u2(0, pi^2) q[0];\nCX q[0], q[1];\nmeasure q -> c;\n"
    )
    assert names(circuit) == ["u2", "cx", "measure", "measure"]
    assert circuit.instructions[0].params == (0.0, math.pi**2)


def test_broadcasting_ranges_sets_and_aliases():
    circuit = parse_qasm(
        source="""
        include "stdgates.inc";
        qubit[4] q;
        qubit a;
        let pair = q[{3, 1}] ++ a;
        cx q[0:1], q[2:3];
        cx a, q[0:2:3];  // start:step:end
        x pair[-1];
        h q[1:-1];
        barrier pair;
        """
    )
    qubits = [tuple(f"{q.register}{q.index}" for q in inst.qubits) for inst in circuit.instructions]
    assert qubits == [
        ("q0", "q2"),
        ("q1", "q3"),
        ("a0", "q0"),
        ("a0", "q2"),
        ("a0",),
        ("q1",),
        ("q2",),
        ("q3",),
        ("q3", "q1", "a0"),
    ]


def test_measure_into_declaration_and_empty_barrier():
    circuit = parse_qasm(source="qubit[2] q; barrier; bit[2] c = measure q;")
    assert circuit.instructions[0] == Instruction("barrier", (Qubit("q", 0), Qubit("q", 1)))
    assert names(circuit) == ["barrier", "measure", "measure"]


def test_includes_relative_to_the_file(tmp_path):
    (tmp_path / "defs.inc").write_text("gate mine a { U(0, 0, 0) a; }\nconst int size = 2;\n")
    path = tmp_path / "main.qasm"
    path.write_text('include "defs.inc";\nqubit[size] q;\nmine q;\n')
    # Gates defined by the program's own files are expanded into their bodies.
    assert names(parse_qasm(path=path)) == ["u", "u"]


def test_gate_definitions_are_expanded():
    circuit = parse_qasm(
        source="""
        OPENQASM 3.0;
        include "stdgates.inc";
        gate twice(t) a, b { rz(2 * t) a; cx a, b; }
        gate outer(t) a, b { h b; twice(t / 2) b, a; }
        qubit[2] q;
        outer(0.6) q[0], q[1];
        """
    )
    assert [(i.name, [q.index for q in i.qubits], i.params) for i in circuit.instructions] == [
        ("h", [1], ()),
        ("rz", [1], (pytest.approx(0.6),)),
        ("cx", [1, 0], ()),
    ]


def test_library_gates_the_compiler_knows_stay_named():
    circuit = parse_qasm(
        source='OPENQASM 2.0; include "qelib1.inc"; qreg q[4]; ccx q[0], q[1], q[2]; c3x q[0], q[1], q[2], q[3];'
    )
    named = names(circuit)
    # ccx is lowered by the compiler itself; c3x comes from its qelib1.inc definition.
    assert named[0] == "ccx" and "c3x" not in named and len(named) > 2


def test_semantic_errors_are_reported_with_locations():
    with pytest.raises(QasmFrontendError, match=r"2:1: undefined gate 'h'"):
        parse_qasm(source="qubit q;\nh q;")
    with pytest.raises(QasmFrontendError, match="expected ';'"):
        parse_qasm(source="qubit q")


@pytest.mark.parametrize(
    ("source", "message"),
    [
        ('include "stdgates.inc"; qubit q; bit b; b = measure q; if (b) x q;', "'if' statement"),
        ('include "stdgates.inc"; qubit[2] q; ctrl @ x q[0], q[1];', "gate modifiers"),
        ('include "stdgates.inc"; x $0;', "hardware qubits"),
        (
            'include "stdgates.inc"; input float theta; qubit q; rx(theta) q;',
            "'theta' is not a compile-time constant",
        ),
        ("int i = 3; i += 1;", "classical variables other than bits"),
        ("qubit q; delay[10ns] q;", "'delay' instruction is not supported"),
        ("int i = 3;", "classical variables other than bits"),
        ("qubit q; measure q;", "must store their result"),
        ("qubit q; ctrl @ gphase(1) q;", "controlled global phases"),
    ],
)
def test_unsupported_features_are_errors_not_silently_dropped(source, message):
    with pytest.raises(QasmFrontendError, match=message):
        parse_qasm(source=source)


@pytest.mark.parametrize("version", ["2.0", "3.0"])
def test_emitted_programs_are_valid_and_round_trip(version):
    source = """
    include "stdgates.inc";
    qubit[3] q;
    bit[3] c;
    U(pi / 2, -3 * pi / 4, 0.125) q[0];
    cx q[0], q[1];
    rz(-pi) q[2];
    barrier q[0], q[2];
    reset q[1];
    c = measure q;
    """
    circuit = parse_qasm(source=source)
    text = emit_qasm(circuit, version=version)
    program = openqasm.parse(text)
    assert program.version == version
    assert openqasm.analyze(program).ok
    assert parse_qasm(source=text).instructions == circuit.instructions


def test_emitted_openqasm2_syntax():
    circuit = parse_qasm(
        source='include "stdgates.inc"; qubit[2] q; bit[2] c; rz(pi / 2) q[0]; c = measure q;'
    )
    assert emit_qasm(circuit) == (
        'OPENQASM 2.0;\ninclude "qelib1.inc";\nqreg q[2];\ncreg c[2];\n'
        "rz(pi / 2) q[0];\nmeasure q[0] -> c[0];\nmeasure q[1] -> c[1];"
    )
    assert "U(pi" in emit_qasm(parse_qasm(source="qubit q; U(pi, 0, 0) q;"), version="3.0")
    angles = parse_qasm(source="qubit q; U(-pi / 2, -3 * pi / 4, 2 * pi) q; U(-0.5, 1e-3, pi) q;")
    assert emit_qasm(angles, version="3.0").splitlines()[-2:] == [
        "U(-pi / 2, -3 * pi / 4, 2 * pi) q[0];",
        "U(-0.5, 0.001, pi) q[0];",
    ]


def test_compile_openqasm3_source_to_native_gates():
    compiled = compile_source('include "stdgates.inc"; qubit[2] q; h q[0]; cx q[0], q[1];')
    assert set(names(compiled)) <= {"rz", "sx", "x", "cx"}


def test_cli_openqasm3_output_and_errors(tmp_path, capsys):
    source = tmp_path / "in.qasm"
    source.write_text('include "stdgates.inc";\nqubit[2] q;\nh q[0];\ncx q[0], q[1];\n')
    out = tmp_path / "out.qasm"
    assert main(["compile", str(source), "-o", str(out), "--qasm-version", "3"]) == 0
    text = out.read_text()
    assert text.startswith('OPENQASM 3.0;\ninclude "stdgates.inc";\nqubit[2] q;')
    assert openqasm.analyze(openqasm.parse(text)).ok
    source.write_text("qubit q;\nnot_a_gate q;\n")
    assert main(["compile", str(source), "-o", str(out)]) == 1
    assert "undefined gate 'not_a_gate'" in capsys.readouterr().err
