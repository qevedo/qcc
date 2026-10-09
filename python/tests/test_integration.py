from qevedo.compiler.compiler import CompileOptions, compile_file
from qevedo.compiler.device import DeviceSpec, default_device
from qevedo.compiler.emit.qasm import emit_qasm
from pathlib import Path


EXAMPLES = Path(__file__).resolve().parents[1] / "examples"
PROFILES = Path(__file__).resolve().parents[1] / "qevedo" / "compiler" / "profiles"


def test_compile_bell_to_native_gates(tmp_path):
    out = tmp_path / "out.qasm"
    circuit = compile_file(EXAMPLES / "bell.qasm", CompileOptions(device=default_device()))
    text = emit_qasm(circuit)
    out.write_text(text)
    assert "h " not in text
    assert "rz(" in text
    assert "cx " in text


def test_cli_stats_device_profile():
    spec = DeviceSpec.from_yaml(PROFILES / "linear_chain.yaml")
    assert spec.topology == "linear_chain"
    assert len(spec.edges) == 15
