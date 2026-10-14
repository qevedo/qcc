//! The hardware model: a modular quantum computer made of several cores.

/// A multi-core quantum computer: `cores` cores of `qubits_per_core` physical
/// qubits each.
///
/// As in Escofet et al. (2023, 2025), every core is all-to-all inside (any
/// two qubits of a core can interact directly) and the cores are connected
/// all-to-all: moving a qubit's state to any other core is one teleportation,
/// one *non-local communication*. A two-qubit gate can only run when both of
/// its qubits are in the same core.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Architecture {
    pub cores: usize,
    pub qubits_per_core: usize,
}

impl Architecture {
    pub fn new(cores: usize, qubits_per_core: usize) -> Result<Architecture, String> {
        if cores < 1 {
            return Err("an architecture needs at least one core".into());
        }
        if qubits_per_core < 2 {
            return Err("each core needs at least two qubits to run a two-qubit gate".into());
        }
        Ok(Architecture {
            cores,
            qubits_per_core,
        })
    }

    /// The total number of physical qubits.
    pub fn capacity(&self) -> usize {
        self.cores * self.qubits_per_core
    }

    /// Check that `qubits` circuit qubits fit.
    pub fn check_fits(&self, qubits: usize) -> Result<(), String> {
        if qubits > self.capacity() {
            return Err(format!(
                "{qubits} qubits do not fit in {} cores of {} qubits",
                self.cores, self.qubits_per_core
            ));
        }
        Ok(())
    }
}
