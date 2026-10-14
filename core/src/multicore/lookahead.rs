//! How strongly a qubit is drawn to a core by the interactions still to come.
//!
//! Escofet et al. weight an interaction `m − t` timeslices ahead by
//! `2^−(m−t)`, so the next slice counts most:
//!
//! ```text
//! w_t(q, q') = Σ_{t < m ≤ T} I(m, q, q') · 2^−(m−t)
//! attr_t(q, c) = Σ_{q' in core c} w_t(q, q')
//! ```
//!
//! A qubit has at most one partner per slice, so `attr_t(q, c)` is a sum over
//! `q`'s future partners. Terms beyond [`HORIZON`] slices are smaller than the
//! rounding error of the sum (at most 1), so they are skipped.

use super::timeslices::Timeslices;

/// Slices looked ahead: 2^−64 is far below double precision relative to 1.
pub const HORIZON: usize = 64;

/// The attraction of qubit `q` to core `core` after slice `t`, given where the
/// other qubits are (`None` for qubits not placed yet).
pub fn attraction(
    slices: &Timeslices,
    t: usize,
    q: usize,
    core: usize,
    core_of: &[Option<usize>],
) -> f64 {
    let mut total = 0.0;
    let mut weight = 1.0;
    for m in t + 1..slices.len().min(t + 1 + HORIZON) {
        weight *= 0.5;
        if let Some(partner) = slices.partner(m, q) {
            if core_of[partner] == Some(core) {
                total += weight;
            }
        }
    }
    total
}

/// The attraction of qubit `q` to every core at once (one pass over its
/// future partners): `out[c] = attraction(slices, t, q, c, core_of)`.
pub fn attraction_to_cores(
    slices: &Timeslices,
    t: usize,
    q: usize,
    core_of: &[Option<usize>],
    out: &mut [f64],
) {
    out.iter_mut().for_each(|x| *x = 0.0);
    let mut weight = 1.0;
    for m in t + 1..slices.len().min(t + 1 + HORIZON) {
        weight *= 0.5;
        if let Some(core) = slices.partner(m, q).and_then(|partner| core_of[partner]) {
            out[core] += weight;
        }
    }
}

/// The look-ahead cost of qubit `q` in every core when the cores are not all
/// linked: `out[c] = Σ w · hops(c, core of the partner)` over its future
/// partners already placed, the teleportations it would need to meet them,
/// weighted like the attraction. With all cores linked to all this is the
/// total weight minus the attraction, so it orders cores the same way.
pub fn future_hops_to_cores(
    slices: &Timeslices,
    t: usize,
    q: usize,
    core_of: &[Option<usize>],
    arch: &super::Architecture,
    out: &mut [f64],
) {
    out.iter_mut().for_each(|x| *x = 0.0);
    let mut weight = 1.0;
    for m in t + 1..slices.len().min(t + 1 + HORIZON) {
        weight *= 0.5;
        if let Some(core) = slices.partner(m, q).and_then(|partner| core_of[partner]) {
            for (c, x) in out.iter_mut().enumerate() {
                *x += weight * arch.hops(c, core) as f64;
            }
        }
    }
}

/// The look-ahead weights `w_t` of every pair of qubits, as a dense matrix.
pub fn weights(slices: &Timeslices, t: usize) -> Vec<Vec<f64>> {
    let n = slices.qubits();
    let mut w = vec![vec![0.0; n]; n];
    let mut weight = 1.0;
    for m in t + 1..slices.len().min(t + 1 + HORIZON) {
        weight *= 0.5;
        for &(a, b) in slices.slice(m) {
            w[a][b] += weight;
            w[b][a] += weight;
        }
    }
    w
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nearer_interactions_weigh_more() {
        // q0 meets q1 one slice after slice 0, q2 two slices after.
        let s = Timeslices::from_pairs(3, &[(0, 1), (0, 1), (0, 2)]).unwrap();
        let core_of = [Some(0), Some(1), Some(2)];
        assert_eq!(attraction(&s, 0, 0, 1, &core_of), 0.5);
        assert_eq!(attraction(&s, 0, 0, 2, &core_of), 0.25);
        assert_eq!(attraction(&s, 0, 0, 0, &core_of), 0.0);
        let mut row = [0.0; 3];
        attraction_to_cores(&s, 0, 0, &core_of, &mut row);
        assert_eq!(row, [0.0, 0.5, 0.25]);
        let w = weights(&s, 0);
        assert_eq!((w[0][1], w[0][2], w[2][0]), (0.5, 0.25, 0.25));
    }
}
