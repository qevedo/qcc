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

/// Map a circuit, given by its two-qubit gates in program order, onto a
/// multi-core architecture with `method` ("hqa", "hqa-no-lookahead",
/// "fgp-roee" or "naive"). Returns (initial assignment, assignment per
/// timeslice, gates per timeslice, non-local communications).
#[pyfunction]
#[allow(clippy::type_complexity)]
fn map_to_cores(
    py: Python<'_>,
    pairs: Vec<(usize, usize)>,
    qubits: usize,
    cores: usize,
    qubits_per_core: usize,
    method: &str,
) -> PyResult<(Vec<usize>, Vec<Vec<usize>>, Vec<Vec<(usize, usize)>>, usize)> {
    use qcc::multicore::{Architecture, Method, Timeslices, map};
    let method = Method::parse(method).map_err(PyValueError::new_err)?;
    let arch = Architecture::new(cores, qubits_per_core).map_err(PyValueError::new_err)?;
    // Slices with more gates than the cores can run at once (odd core sizes) are split.
    let slices = Timeslices::from_pairs(qubits, &pairs)
        .map_err(PyValueError::new_err)?
        .fitted_to(&arch);
    let mapping = py
        .detach(|| map(&slices, &arch, method))
        .map_err(PyValueError::new_err)?;
    let gates = (0..slices.len())
        .map(|t| slices.slice(t).to_vec())
        .collect();
    let communications = mapping.communications();
    Ok((mapping.initial, mapping.slices, gates, communications))
}

/// Lower and upper bounds on non-local communications for a random circuit
/// with `two_qubit_gates` gates on `qubits` qubits over `cores` cores.
#[pyfunction]
fn random_circuit_bounds(cores: usize, qubits: usize, two_qubit_gates: usize) -> (f64, f64) {
    qcc::multicore::random_circuit_bounds(cores, qubits, two_qubit_gates)
}

#[pymodule]
fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(compile, m)?)?;
    m.add_function(wrap_pyfunction!(run_passes, m)?)?;
    m.add_function(wrap_pyfunction!(map_to_cores, m)?)?;
    m.add_function(wrap_pyfunction!(random_circuit_bounds, m)?)?;
    Ok(())
}
