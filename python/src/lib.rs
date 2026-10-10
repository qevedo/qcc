//! Python bindings of the QCC core, imported as `qevedo.compiler._core`.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

/// An instruction as Python passes it: (name, qubits, params, clbits).
type Raw = (String, Vec<usize>, Vec<f64>, Vec<usize>);

/// Lower `instructions` to the `native` gates and optimize them with the
/// default pipeline. Raises ValueError when a gate cannot be lowered.
#[pyfunction]
fn compile(py: Python<'_>, instructions: Vec<Raw>, native: Vec<String>) -> PyResult<Vec<Raw>> {
    let instructions: Vec<qcc::Instruction> = instructions
        .into_iter()
        .map(|(name, qubits, params, clbits)| qcc::Instruction {
            name: name.to_lowercase(),
            qubits,
            params,
            clbits,
        })
        .collect();
    let out = py
        .detach(move || qcc::compile(instructions, native))
        .map_err(PyValueError::new_err)?;
    Ok(out
        .into_iter()
        .map(|i| (i.name, i.qubits, i.params, i.clbits))
        .collect())
}

/// Run the named passes of the core in order ("decompose", "merge_1q",
/// "commutative_cancellation", "resynthesize_2q"), for testing.
#[pyfunction]
fn run_passes(
    py: Python<'_>,
    instructions: Vec<Raw>,
    native: Vec<String>,
    passes: Vec<String>,
) -> PyResult<Vec<Raw>> {
    let instructions: Vec<qcc::Instruction> = instructions
        .into_iter()
        .map(|(name, qubits, params, clbits)| qcc::Instruction {
            name: name.to_lowercase(),
            qubits,
            params,
            clbits,
        })
        .collect();
    let out = py
        .detach(move || {
            let passes: Vec<&str> = passes.iter().map(String::as_str).collect();
            qcc::passes::run_passes(instructions, native, &passes)
        })
        .map_err(PyValueError::new_err)?;
    Ok(out
        .into_iter()
        .map(|i| (i.name, i.qubits, i.params, i.clbits))
        .collect())
}

#[pymodule]
fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(compile, m)?)?;
    m.add_function(wrap_pyfunction!(run_passes, m)?)?;
    Ok(())
}
