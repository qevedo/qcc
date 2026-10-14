//! Map a random circuit onto a multi-core quantum computer with each method
//! and compare the non-local communications they need.
//!
//!     cargo run --release --example multicore -- [qubits] [cores]
//!
//! The circuit has 2 × qubits layers, each pairing up all qubits at random,
//! like the random circuits of Escofet et al.

use qcc::multicore::{Architecture, Method, Timeslices, map, random_circuit_bounds};
use qcc::rng::Rng;
use std::time::Instant;

fn main() {
    let args: Vec<usize> = std::env::args()
        .skip(1)
        .filter_map(|a| a.parse().ok())
        .collect();
    let qubits = args.first().copied().unwrap_or(100);
    let cores = args.get(1).copied().unwrap_or(10);
    let arch = Architecture::new(cores, qubits.div_ceil(cores)).expect("a valid architecture");

    let mut rng = Rng::new(1);
    let mut pairs = Vec::new();
    for _ in 0..2 * qubits {
        let mut order: Vec<usize> = (0..qubits).collect();
        for i in (1..qubits).rev() {
            order.swap(i, rng.below(i + 1));
        }
        pairs.extend(
            order
                .chunks(2)
                .filter(|c| c.len() == 2)
                .map(|c| (c[0], c[1])),
        );
    }
    let slices = Timeslices::from_pairs(qubits, &pairs)
        .expect("valid gates")
        .fitted_to(&arch);
    let (lower, upper) = random_circuit_bounds(cores, qubits, pairs.len());
    println!(
        "{qubits} qubits on {cores} cores of {}: {} two-qubit gates in {} timeslices",
        arch.qubits_per_core,
        pairs.len(),
        slices.len()
    );
    println!("random-circuit bounds: {lower:.0} to {upper:.0} communications\n");
    println!("{:<20} {:>15} {:>10}", "method", "communications", "time");
    for method in [
        Method::Hqa,
        Method::HqaWithoutLookahead,
        Method::FgpRoee,
        Method::Naive,
    ] {
        let start = Instant::now();
        let mapping = map(&slices, &arch, method).expect("a mapping");
        let elapsed = start.elapsed();
        mapping.validate(&arch, &slices).expect("a valid mapping");
        println!(
            "{:<20} {:>15} {:>10.1?}",
            format!("{method:?}"),
            mapping.communications(),
            elapsed
        );
    }
}
