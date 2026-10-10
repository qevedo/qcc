//! The instruction list the passes work on: gates, measurements, resets and
//! barriers on numbered qubits and classical bits.

#[derive(Clone, Debug, PartialEq)]
pub struct Instruction {
    pub name: String,
    pub qubits: Vec<usize>,
    pub params: Vec<f64>,
    pub clbits: Vec<usize>,
}

/// Operations that are not gates: they pass through unchanged and block optimizations.
pub const NON_GATES: [&str; 3] = ["barrier", "measure", "reset"];

impl Instruction {
    pub fn gate(name: &str, qubits: Vec<usize>, params: Vec<f64>) -> Instruction {
        Instruction {
            name: name.to_string(),
            qubits,
            params,
            clbits: vec![],
        }
    }

    pub fn is_gate(&self) -> bool {
        !NON_GATES.contains(&self.name.as_str()) && self.clbits.is_empty()
    }

    /// The same instruction with qubit `i` replaced by `qubits[i]`.
    pub fn relabel(&self, qubits: &[usize]) -> Instruction {
        Instruction {
            name: self.name.clone(),
            qubits: self.qubits.iter().map(|&q| qubits[q]).collect(),
            params: self.params.clone(),
            clbits: self.clbits.clone(),
        }
    }
}

/// The cost the passes minimise: two-qubit gates first, then all gates.
pub fn cost(instructions: &[Instruction]) -> (usize, usize) {
    (
        instructions.iter().filter(|i| i.qubits.len() > 1).count(),
        instructions.len(),
    )
}
