//! Assignments of qubits to cores, and what moving between them costs.

use super::architecture::Architecture;
use super::timeslices::Timeslices;

/// Which core each circuit qubit is in: `core_of[q]`.
pub type Assignment = Vec<usize>;

/// The result of a mapper: where every qubit is before the circuit starts,
/// and where it is during each timeslice.
#[derive(Clone, Debug)]
pub struct Mapping {
    pub initial: Assignment,
    /// `slices[t]`: the assignment while slice `t` runs.
    pub slices: Vec<Assignment>,
}

impl Mapping {
    /// The number of non-local communications: every time a qubit is in a
    /// different core than in the previous step, its state was teleported
    /// once. The initial placement is free.
    pub fn communications(&self) -> usize {
        let mut previous = &self.initial;
        let mut total = 0;
        for next in &self.slices {
            total += moves(previous, next);
            previous = next;
        }
        total
    }

    /// Check that every assignment respects the core capacities and that every
    /// gate of each slice has both qubits in the same core.
    pub fn validate(&self, arch: &Architecture, slices: &Timeslices) -> Result<(), String> {
        if self.slices.len() != slices.len() {
            return Err(format!(
                "{} assignments for {} timeslices",
                self.slices.len(),
                slices.len()
            ));
        }
        check_capacity(arch, &self.initial).map_err(|e| format!("initial assignment: {e}"))?;
        for (t, assignment) in self.slices.iter().enumerate() {
            check_capacity(arch, assignment).map_err(|e| format!("timeslice {t}: {e}"))?;
            for &(a, b) in slices.slice(t) {
                if assignment[a] != assignment[b] {
                    return Err(format!(
                        "timeslice {t}: qubits {a} and {b} interact but are in cores {} and {}",
                        assignment[a], assignment[b]
                    ));
                }
            }
        }
        Ok(())
    }
}

/// How many qubits are in a different core in `next` than in `previous`.
pub fn moves(previous: &[usize], next: &[usize]) -> usize {
    previous.iter().zip(next).filter(|(a, b)| a != b).count()
}

/// How many qubits each core holds.
pub fn occupancy(arch: &Architecture, assignment: &[usize]) -> Vec<usize> {
    let mut count = vec![0; arch.cores];
    for &core in assignment {
        count[core] += 1;
    }
    count
}

fn check_capacity(arch: &Architecture, assignment: &[usize]) -> Result<(), String> {
    for (core, &count) in occupancy(arch, assignment).iter().enumerate() {
        if count > arch.qubits_per_core {
            return Err(format!(
                "core {core} holds {count} qubits, more than its {}",
                arch.qubits_per_core
            ));
        }
    }
    Ok(())
}

/// Fill the cores in order: qubits `0..k` in core 0, the next `k` in core 1, …
pub fn in_order(arch: &Architecture, qubits: usize) -> Assignment {
    (0..qubits).map(|q| q / arch.qubits_per_core).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn communications_count_moved_qubits() {
        let mapping = Mapping {
            initial: vec![0, 0, 1, 1],
            slices: vec![vec![0, 1, 0, 1], vec![0, 1, 0, 1], vec![1, 1, 0, 0]],
        };
        // Two moves, then none, then two.
        assert_eq!(mapping.communications(), 4);
    }

    #[test]
    fn validation_catches_split_gates_and_full_cores() {
        let arch = Architecture::new(2, 2).unwrap();
        let slices = Timeslices::from_pairs(4, &[(0, 1)]).unwrap();
        let ok = Mapping {
            initial: vec![0, 0, 1, 1],
            slices: vec![vec![0, 0, 1, 1]],
        };
        assert!(ok.validate(&arch, &slices).is_ok());
        let split = Mapping {
            initial: vec![0, 0, 1, 1],
            slices: vec![vec![0, 1, 0, 1]],
        };
        assert!(
            split
                .validate(&arch, &slices)
                .unwrap_err()
                .contains("interact")
        );
        let full = Mapping {
            initial: vec![0, 0, 1, 1],
            slices: vec![vec![0, 0, 0, 1]],
        };
        assert!(
            full.validate(&arch, &slices)
                .unwrap_err()
                .contains("more than")
        );
    }
}
