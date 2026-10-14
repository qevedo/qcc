//! A circuit as the mappers see it: its two-qubit gates, grouped into
//! timeslices.

use crate::ir::Instruction;

/// The two-qubit interactions of a circuit, grouped into timeslices.
///
/// A timeslice is a set of two-qubit gates that can run at the same time, so
/// each qubit takes part in at most one gate per slice. Slices are built as
/// soon as possible: each gate goes into the first slice after the last slice
/// of either of its qubits. Single-qubit gates do not affect which core a
/// qubit must be in, so they are left out.
#[derive(Clone, Debug)]
pub struct Timeslices {
    qubits: usize,
    slices: Vec<Vec<(usize, usize)>>,
    /// `partners[t][q]`: the qubit that `q` interacts with in slice `t`.
    partners: Vec<Vec<Option<usize>>>,
}

impl Timeslices {
    /// Slices from the two-qubit gates of a circuit, in program order.
    pub fn from_pairs(qubits: usize, gates: &[(usize, usize)]) -> Result<Timeslices, String> {
        let mut last = vec![0usize; qubits];
        let mut slices: Vec<Vec<(usize, usize)>> = Vec::new();
        for &(a, b) in gates {
            if a >= qubits || b >= qubits {
                return Err(format!(
                    "gate on ({a}, {b}) uses a qubit beyond the {qubits} of the circuit"
                ));
            }
            if a == b {
                return Err(format!("gate on ({a}, {b}) uses the same qubit twice"));
            }
            let t = last[a].max(last[b]);
            if t == slices.len() {
                slices.push(Vec::new());
            }
            slices[t].push((a, b));
            last[a] = t + 1;
            last[b] = t + 1;
        }
        Ok(Timeslices::from_slices(qubits, slices))
    }

    /// Slices from a compiled instruction list on qubits `0..qubits`.
    ///
    /// Gates on three or more qubits must be decomposed first; barriers,
    /// measurements and resets do not constrain the mapping.
    pub fn from_instructions(
        qubits: usize,
        instructions: &[Instruction],
    ) -> Result<Timeslices, String> {
        let mut pairs = Vec::new();
        for inst in instructions {
            if !inst.is_gate() {
                continue;
            }
            match inst.qubits.len() {
                0 | 1 => {}
                2 => pairs.push((inst.qubits[0], inst.qubits[1])),
                n => {
                    return Err(format!(
                        "'{}' acts on {n} qubits: decompose the circuit into one- and two-qubit gates first",
                        inst.name
                    ));
                }
            }
        }
        Timeslices::from_pairs(qubits, &pairs)
    }

    /// These slices with any slice that has more gates than `arch` can run at
    /// once (each core holds ⌊size / 2⌋ gates) split into consecutive slices.
    /// Only cores of odd size ever need it.
    pub fn fitted_to(&self, arch: &super::Architecture) -> Timeslices {
        let most = arch.cores * (arch.qubits_per_core / 2);
        if self.slices.iter().all(|s| s.len() <= most) {
            return self.clone();
        }
        let chunks: Vec<Vec<(usize, usize)>> = self
            .slices
            .iter()
            .flat_map(|s| s.chunks(most).map(<[_]>::to_vec))
            .collect();
        Timeslices::from_slices(self.qubits, chunks)
    }

    fn from_slices(qubits: usize, slices: Vec<Vec<(usize, usize)>>) -> Timeslices {
        let partners = slices
            .iter()
            .map(|slice| {
                let mut partner = vec![None; qubits];
                for &(a, b) in slice {
                    partner[a] = Some(b);
                    partner[b] = Some(a);
                }
                partner
            })
            .collect();
        Timeslices {
            qubits,
            slices,
            partners,
        }
    }

    /// Whether every slice can run at once on `arch`.
    pub fn fits(&self, arch: &super::Architecture) -> bool {
        let most = arch.cores * (arch.qubits_per_core / 2);
        self.slices.iter().all(|s| s.len() <= most)
    }

    pub fn qubits(&self) -> usize {
        self.qubits
    }

    pub fn len(&self) -> usize {
        self.slices.len()
    }

    pub fn is_empty(&self) -> bool {
        self.slices.is_empty()
    }

    /// The gates of slice `t`.
    pub fn slice(&self, t: usize) -> &[(usize, usize)] {
        &self.slices[t]
    }

    /// The qubit that `q` interacts with in slice `t`, if any.
    pub fn partner(&self, t: usize, q: usize) -> Option<usize> {
        self.partners[t][q]
    }

    /// The total number of two-qubit gates.
    pub fn gates(&self) -> usize {
        self.slices.iter().map(Vec::len).sum()
    }

    /// How many times each pair of qubits interacts over the whole circuit, as
    /// a dense symmetric matrix: the "total interaction graph".
    pub fn interaction_counts(&self) -> Vec<Vec<f64>> {
        let mut w = vec![vec![0.0; self.qubits]; self.qubits];
        for slice in &self.slices {
            for &(a, b) in slice {
                w[a][b] += 1.0;
                w[b][a] += 1.0;
            }
        }
        w
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gates_go_into_the_earliest_slice() {
        // q0-q1 and q2-q3 run together; q1-q2 waits for both.
        let s = Timeslices::from_pairs(4, &[(0, 1), (2, 3), (1, 2), (0, 3)]).unwrap();
        assert_eq!(s.len(), 2);
        assert_eq!(s.slice(0), &[(0, 1), (2, 3)]);
        assert_eq!(s.slice(1), &[(1, 2), (0, 3)]);
        assert_eq!(s.partner(1, 3), Some(0));
        assert_eq!(s.partner(0, 2), Some(3));
        assert_eq!(s.gates(), 4);
    }

    #[test]
    fn slices_too_wide_for_odd_cores_are_split() {
        // Three simultaneous gates, two cores of three qubits: at most two at once.
        let s = Timeslices::from_pairs(6, &[(0, 3), (1, 4), (2, 5)]).unwrap();
        let arch = super::super::Architecture::new(2, 3).unwrap();
        assert!(!s.fits(&arch));
        let fitted = s.fitted_to(&arch);
        assert!(fitted.fits(&arch));
        assert_eq!((fitted.len(), fitted.gates()), (2, 3));
    }

    #[test]
    fn instructions_skip_single_qubit_gates_and_reject_wide_ones() {
        let g = |name: &str, q: &[usize]| Instruction::gate(name, q.to_vec(), vec![]);
        let s = Timeslices::from_instructions(
            3,
            &[
                g("h", &[0]),
                g("cx", &[0, 1]),
                g("x", &[2]),
                g("cz", &[1, 2]),
            ],
        )
        .unwrap();
        assert_eq!(s.len(), 2);
        assert!(
            Timeslices::from_instructions(3, &[g("ccx", &[0, 1, 2])])
                .unwrap_err()
                .contains("decompose")
        );
        assert!(Timeslices::from_pairs(2, &[(0, 2)]).is_err());
    }
}
