//! Routing: making every two-qubit gate act on qubits the chip couples.
//!
//! A chip couples only some pairs of its qubits ([`CouplingGraph`]: a line, a
//! ring, a grid, or any connected set of couplings). A router picks where
//! each logical qubit starts (its *layout*) and inserts SWAP gates so that
//! every two-qubit gate runs on coupled qubits, with as few SWAPs as it can.
//!
//! [`layout_and_route`] chooses the layout and routes with SABRE
//! ([`sabre`]); [`route`] routes from a given layout, and [`route_ahead`]
//! also prepares for gates that come after the circuit. Both return a
//! [`Routed`] circuit on physical qubits with the initial and final layouts.
//!
//! | file | what it holds |
//! |---|---|
//! | `coupling.rs` | coupling graphs and their distances |
//! | `sabre.rs` | the SABRE router and its layout search |
//! | `perfect.rs` | layouts that need no SWAP, when the circuit's pairs fit the chip |

mod coupling;
pub mod perfect;
pub mod sabre;

pub use coupling::CouplingGraph;
pub use sabre::{Routed, SabreOptions, layout_and_route, route, route_ahead};
