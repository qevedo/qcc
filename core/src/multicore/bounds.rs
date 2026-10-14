//! Theoretical bounds on non-local communications for random circuits, from
//! Escofet et al. (ACM TQC 2025, Theorems 3.2 and 3.3).
//!
//! For a random circuit on `q` qubits spread over `N` cores, with `g` gates of
//! which a fraction `f` act on two qubits, a mapper that keeps qubits in place
//! needs about `(N − 1) g f q / (N (q − 1))` communications at least (each
//! two-qubit gate lands across cores with probability `(N − 1) q / (N (q − 1))`
//! and costs at least one teleportation), and the naive strategy about twice
//! that.

/// `(lower, upper)` bounds for `two_qubit_gates` (= g·f) random gates on
/// `qubits` qubits over `cores` cores.
pub fn random_circuit_bounds(cores: usize, qubits: usize, two_qubit_gates: usize) -> (f64, f64) {
    let (n, q, gf) = (cores as f64, qubits as f64, two_qubit_gates as f64);
    let lower = (n - 1.0) * gf * q / (n * (q - 1.0));
    (lower, 2.0 * lower)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn one_core_needs_no_communication() {
        assert_eq!(random_circuit_bounds(1, 10, 100), (0.0, 0.0));
        let (lower, upper) = random_circuit_bounds(10, 100, 1000);
        assert!((lower - 909.0909).abs() < 1e-3 && (upper - 2.0 * lower).abs() < 1e-9);
    }
}
