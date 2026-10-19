//! QCC core: gate synthesis and optimization passes of the Quantum Compiler
//! Collection, used by the `qevedo-compiler` Python package.

// Index loops read best in the small fixed-size matrix code.
#![allow(clippy::needless_range_loop)]

pub mod decompose;
pub mod device;
pub mod gates;
pub mod ir;
pub mod linalg;
pub mod lower;
pub mod multicore;
pub mod numeric;
pub mod one_qubit;
pub mod passes;
pub mod permutation;
pub mod phase_folding;
pub mod rng;
pub mod routing;
pub mod symmetric_blocks;
#[cfg(test)]
pub(crate) mod testing;
pub mod two_qubit;

pub use ir::Instruction;
pub use passes::compile;
