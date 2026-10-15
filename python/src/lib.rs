//! Python bindings of the QCC core, imported as `qevedo.compiler._core`.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

/// An instruction as Python passes it: (name, qubits, params, clbits).
type Raw = (String, Vec<usize>, Vec<f64>, Vec<usize>);

fn from_raw(instructions: Vec<Raw>) -> Vec<qcc::Instruction> {
    instructions
        .into_iter()
        .map(|(name, qubits, params, clbits)| qcc::Instruction {
            name: name.to_lowercase(),
            qubits,
            params,
            clbits,
        })
        .collect()
}

/// A timing model as Python passes it: (one_qubit, two_qubit, swap, measure,
/// epr_pair, teleport, link_capacity).
type RawTiming = (f64, f64, f64, f64, f64, f64, usize);

fn timing(raw: Option<RawTiming>) -> qcc::multicore::Timing {
    match raw {
        None => qcc::multicore::Timing::default(),
        Some((one_qubit, two_qubit, swap, measure, epr_pair, teleport, link_capacity)) => {
            qcc::multicore::Timing {
                one_qubit,
                two_qubit,
                swap,
                measure,
                epr_pair,
                teleport,
                link_capacity,
            }
        }
    }
}

fn to_raw(instructions: Vec<qcc::Instruction>) -> Vec<Raw> {
    instructions
        .into_iter()
        .map(|i| (i.name, i.qubits, i.params, i.clbits))
        .collect()
}

/// Lower `instructions` to the `native` gates and optimize them with the
/// default pipeline. Raises ValueError when a gate cannot be lowered.
#[pyfunction]
fn compile(py: Python<'_>, instructions: Vec<Raw>, native: Vec<String>) -> PyResult<Vec<Raw>> {
    let instructions = from_raw(instructions);
    let out = py
        .detach(move || qcc::compile(instructions, native))
        .map_err(PyValueError::new_err)?;
    Ok(to_raw(out))
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
    let instructions = from_raw(instructions);
    let out = py
        .detach(move || {
            let passes: Vec<&str> = passes.iter().map(String::as_str).collect();
            qcc::passes::run_passes(instructions, native, &passes)
        })
        .map_err(PyValueError::new_err)?;
    Ok(to_raw(out))
}

/// Map a circuit, given by its two-qubit gates in program order, onto a
/// multi-core architecture with `method` ("hqa", "hqa-no-lookahead",
/// "hqa-uniform-cost", "fgp-roee" or "naive"). The cores are linked by `links`
/// (pairs of cores), or all to all when it is None; "hqa-latency" plans for
/// `timing` (the default timing when None). Returns (initial assignment,
/// assignment per timeslice, gates per timeslice, communications, moved
/// qubits).
#[pyfunction]
#[pyo3(signature = (pairs, qubits, cores, qubits_per_core, method, links=None, timing=None))]
#[allow(clippy::type_complexity, clippy::too_many_arguments)]
fn map_to_cores(
    py: Python<'_>,
    pairs: Vec<(usize, usize)>,
    qubits: usize,
    cores: usize,
    qubits_per_core: usize,
    method: &str,
    links: Option<Vec<(usize, usize)>>,
    timing: Option<RawTiming>,
) -> PyResult<(
    Vec<usize>,
    Vec<Vec<usize>>,
    Vec<Vec<(usize, usize)>>,
    usize,
    usize,
)> {
    use qcc::multicore::{Architecture, Method, Timeslices, map_timed};
    let timing = self::timing(timing);
    let method = Method::parse(method).map_err(PyValueError::new_err)?;
    let arch = match links {
        Some(links) => Architecture::with_links(cores, qubits_per_core, &links),
        None => Architecture::new(cores, qubits_per_core),
    }
    .map_err(PyValueError::new_err)?;
    // Slices with more gates than the cores can run at once (odd core sizes) are split.
    let slices = Timeslices::from_pairs(qubits, &pairs)
        .map_err(PyValueError::new_err)?
        .fitted_to(&arch);
    let mapping = py
        .detach(|| map_timed(&slices, &arch, method, &timing))
        .map_err(PyValueError::new_err)?;
    let gates = (0..slices.len())
        .map(|t| slices.slice(t).to_vec())
        .collect();
    let (communications, moved) = (mapping.communications(), mapping.moved_qubits());
    Ok((
        mapping.initial,
        mapping.slices,
        gates,
        communications,
        moved,
    ))
}

/// Route a circuit (instructions on logical qubits `0..logical`) onto a chip of
/// `physical` qubits coupled by `edges`, with SABRE (layout search included).
/// Returns (instructions on physical qubits, initial layout, final layout,
/// SWAPs inserted).
#[pyfunction]
#[pyo3(signature = (instructions, logical, physical, edges, trials=8, seed=1))]
#[allow(clippy::type_complexity)]
fn route(
    py: Python<'_>,
    instructions: Vec<Raw>,
    logical: usize,
    physical: usize,
    edges: Vec<(usize, usize)>,
    trials: usize,
    seed: u64,
) -> PyResult<(Vec<Raw>, Vec<usize>, Vec<usize>, usize)> {
    use qcc::routing::{CouplingGraph, SabreOptions, layout_and_route};
    let graph = CouplingGraph::new(physical, &edges).map_err(PyValueError::new_err)?;
    let instructions = from_raw(instructions);
    let options = SabreOptions {
        trials,
        seed,
        ..SabreOptions::default()
    };
    let routed = py
        .detach(|| layout_and_route(&instructions, logical, &graph, &options))
        .map_err(PyValueError::new_err)?;
    Ok((
        to_raw(routed.instructions),
        routed.initial_layout,
        routed.final_layout,
        routed.swaps,
    ))
}

/// Map a circuit onto cores with `method`, as `map_to_cores`, then route it
/// inside each core, whose qubits are coupled by `chip_edges`, with SABRE.
/// Returns (initial assignment, assignment per timeslice, gates per
/// timeslice, communications, moved qubits, instructions on physical qubits,
/// initial layout, final layout, SWAPs, teleportations).
#[pyfunction]
#[pyo3(signature = (instructions, qubits, cores, qubits_per_core, method, links, chip_edges, timing=None, trials=8, seed=1, partner_slots=true, look_ahead=4))]
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
fn map_and_route(
    py: Python<'_>,
    instructions: Vec<Raw>,
    qubits: usize,
    cores: usize,
    qubits_per_core: usize,
    method: &str,
    links: Vec<(usize, usize)>,
    chip_edges: Vec<(usize, usize)>,
    timing: Option<RawTiming>,
    trials: usize,
    seed: u64,
    partner_slots: bool,
    look_ahead: usize,
) -> PyResult<(
    Vec<usize>,
    Vec<Vec<usize>>,
    Vec<Vec<(usize, usize)>>,
    usize,
    usize,
    Vec<Raw>,
    Vec<usize>,
    Vec<usize>,
    usize,
    usize,
)> {
    use qcc::multicore::{Architecture, Method, Timeslices};
    use qcc::routing::{CouplingGraph, SabreOptions};
    let method = Method::parse(method).map_err(PyValueError::new_err)?;
    let arch =
        Architecture::with_links(cores, qubits_per_core, &links).map_err(PyValueError::new_err)?;
    let chip = CouplingGraph::new(qubits_per_core, &chip_edges).map_err(PyValueError::new_err)?;
    let instructions = from_raw(instructions);
    let options = qcc::multicore::PhysicalOptions {
        sabre: SabreOptions {
            trials,
            seed,
            ..SabreOptions::default()
        },
        partner_slots,
        look_ahead,
    };
    let (mapping, physical) = py
        .detach(|| {
            qcc::multicore::map_and_route(
                &instructions,
                qubits,
                &arch,
                &chip,
                method,
                &self::timing(timing),
                &options,
            )
        })
        .map_err(PyValueError::new_err)?;
    let slices = Timeslices::from_instructions(qubits, &instructions)
        .map_err(PyValueError::new_err)?
        .fitted_to(&arch);
    let gates = (0..slices.len())
        .map(|t| slices.slice(t).to_vec())
        .collect();
    let (communications, moved) = (mapping.communications(), mapping.moved_qubits());
    Ok((
        mapping.initial,
        mapping.slices,
        gates,
        communications,
        moved,
        to_raw(physical.instructions),
        physical.initial_layout,
        physical.final_layout,
        physical.swaps,
        physical.teleportations,
    ))
}

/// Schedule a multi-core circuit on physical qubits (as `map_and_route`
/// returns it) on cores linked by `links`, under `timing` (see
/// `qcc::multicore::Timing`). Returns (start of each instruction, duration,
/// time teleportations waited for EPR pairs, hops per link).
#[pyfunction]
#[allow(clippy::type_complexity)]
fn schedule(
    instructions: Vec<Raw>,
    cores: usize,
    qubits_per_core: usize,
    links: Vec<(usize, usize)>,
    timing: RawTiming,
) -> PyResult<(Vec<f64>, f64, f64, Vec<((usize, usize), usize)>)> {
    let arch = qcc::multicore::Architecture::with_links(cores, qubits_per_core, &links)
        .map_err(PyValueError::new_err)?;
    let timing = self::timing(Some(timing));
    let s = qcc::multicore::schedule(&from_raw(instructions), &arch, &timing)
        .map_err(PyValueError::new_err)?;
    Ok((s.start, s.duration, s.epr_wait, s.link_hops))
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
    m.add_function(wrap_pyfunction!(route, m)?)?;
    m.add_function(wrap_pyfunction!(map_and_route, m)?)?;
    m.add_function(wrap_pyfunction!(schedule, m)?)?;
    Ok(())
}
