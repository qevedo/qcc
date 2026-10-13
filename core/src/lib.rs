//! QCC core: gate synthesis and optimization passes of the Quantum Compiler
//! Collection, used by the `qevedo-compiler` Python package.

// Index loops read best in the small fixed-size matrix code.
#![allow(clippy::needless_range_loop)]

pub mod decompose;
pub mod gates;
pub mod ir;
pub mod linalg;
pub mod lower;
pub mod numeric;
pub mod one_qubit;
pub mod passes;
pub mod rng;
#[cfg(test)]
pub(crate) mod testing;
pub mod two_qubit;

pub use ir::Instruction;
pub use passes::compile;
