//! Overall Extreme Exchange (OEE): balanced k-way graph partitioning.
//!
//! An extension of Kernighan–Lin to k parts of fixed sizes: in
//! each pass it swaps, one pair at a time, the two nodes in different parts
//! whose exchange most reduces the weight of edges between parts, locking
//! swapped nodes, and keeps the best prefix of that sequence of swaps. Passes
//! repeat until a pass gains nothing.
//!
//! The mappers use it twice: to place qubits before the circuit starts (on
//! the total interaction graph), and, in relaxed form (rOEE, stop as soon as
//! the partition is valid), inside the FGP-rOEE baseline.

/// How many best candidates per side are tried when looking for the best
/// swap between two parts. Exact for the best pair unless the edge between
/// the two candidates reorders them, which is rare and only affects quality.
const CANDIDATES: usize = 8;

/// A test that a partition is acceptable, for relaxed OEE.
pub type Validity<'a> = &'a dyn Fn(&[usize]) -> bool;

/// Partition `nodes` nodes into the parts given by `part` (sizes are kept)
/// to minimise the total weight `w` of edges between parts. `w` is a dense
/// symmetric matrix over at least `nodes` nodes. `valid`, when given, ends the
/// search as soon as it accepts the current partition (relaxed OEE).
pub fn partition(w: &[Vec<f64>], part: &mut [usize], parts: usize, valid: Option<Validity>) {
    let n = part.len();
    if valid.is_some_and(|ok| ok(part)) {
        return;
    }
    // d[v][p]: the weight between node v and part p.
    let mut d = vec![vec![0.0; parts]; n];
    for passes in 0.. {
        if passes > 4 * n {
            break;
        }
        for (v, row) in d.iter_mut().enumerate() {
            row.iter_mut().for_each(|x| *x = 0.0);
            for u in 0..n {
                if u != v {
                    row[part[u]] += w[v][u];
                }
            }
        }
        let mut locked = vec![false; n];
        let mut swaps: Vec<(usize, usize)> = Vec::new();
        let mut gains: Vec<f64> = Vec::new();
        let mut became_valid = None;
        while let Some((a, b, gain)) = best_swap(w, part, &d, &locked, parts) {
            apply_swap(w, part, &mut d, a, b);
            locked[a] = true;
            locked[b] = true;
            swaps.push((a, b));
            gains.push(gain);
            if valid.is_some_and(|ok| ok(part)) {
                became_valid = Some(swaps.len());
                break;
            }
        }
        // Keep the prefix with the largest total gain (or, in relaxed mode,
        // the first prefix that made the partition valid).
        let keep = match became_valid {
            Some(k) => k,
            None => {
                let (mut best, mut best_k, mut running) = (0.0, 0, 0.0);
                for (k, g) in gains.iter().enumerate() {
                    running += g;
                    if running > best + 1e-12 {
                        (best, best_k) = (running, k + 1);
                    }
                }
                best_k
            }
        };
        for &(a, b) in swaps[keep..].iter().rev() {
            part.swap(a, b);
        }
        if became_valid.is_some() || keep == 0 {
            break;
        }
    }
}

/// The unlocked pair in different parts whose swap gains the most.
fn best_swap(
    w: &[Vec<f64>],
    part: &[usize],
    d: &[Vec<f64>],
    locked: &[bool],
    parts: usize,
) -> Option<(usize, usize, f64)> {
    // The unlocked nodes of each part.
    let mut members = vec![Vec::new(); parts];
    for (v, &p) in part.iter().enumerate() {
        if !locked[v] {
            members[p].push(v);
        }
    }
    let mut best: Option<(usize, usize, f64)> = None;
    for pa in 0..parts {
        for pb in pa + 1..parts {
            // Single-node gains of moving into the other part, best first.
            let top = |from: usize, to: usize| {
                let mut c: Vec<(f64, usize)> = members[from]
                    .iter()
                    .copied()
                    .map(|v| (d[v][to] - d[v][from], v))
                    .collect();
                let order =
                    |x: &(f64, usize), y: &(f64, usize)| y.0.total_cmp(&x.0).then(x.1.cmp(&y.1));
                if c.len() > CANDIDATES {
                    c.select_nth_unstable_by(CANDIDATES - 1, order);
                    c.truncate(CANDIDATES);
                }
                c.sort_by(order);
                c
            };
            for &(ga, a) in &top(pa, pb) {
                for &(gb, b) in &top(pb, pa) {
                    let gain = ga + gb - 2.0 * w[a][b];
                    if best.is_none_or(|(_, _, g)| gain > g + 1e-12) {
                        best = Some((a, b, gain));
                    }
                }
            }
        }
    }
    best
}

fn apply_swap(w: &[Vec<f64>], part: &mut [usize], d: &mut [Vec<f64>], a: usize, b: usize) {
    let (pa, pb) = (part[a], part[b]);
    for (v, row) in d.iter_mut().enumerate() {
        if v != a {
            row[pa] -= w[v][a];
            row[pb] += w[v][a];
        }
        if v != b {
            row[pb] -= w[v][b];
            row[pa] += w[v][b];
        }
    }
    part.swap(a, b);
}

/// The total weight of edges between different parts.
pub fn cut(w: &[Vec<f64>], part: &[usize]) -> f64 {
    let n = part.len();
    (0..n)
        .flat_map(|a| (a + 1..n).map(move |b| (a, b)))
        .filter(|&(a, b)| part[a] != part[b])
        .map(|(a, b)| w[a][b])
        .sum()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn finds_the_obvious_clusters() {
        // Two cliques of four, interleaved in the starting partition.
        let n = 8;
        let mut w = vec![vec![0.0; n]; n];
        for a in 0..n {
            for b in 0..n {
                if a != b && a % 2 == b % 2 {
                    w[a][b] = 1.0;
                }
            }
        }
        w[0][1] = 0.1;
        w[1][0] = 0.1;
        let mut part = vec![0, 0, 0, 0, 1, 1, 1, 1];
        partition(&w, &mut part, 2, None);
        assert_eq!(cut(&w, &part), 0.1);
        assert_eq!(
            part.iter().filter(|&&p| p == 0).count(),
            4,
            "sizes are kept"
        );
    }

    #[test]
    fn relaxed_mode_stops_when_valid() {
        let n = 4;
        let mut w = vec![vec![0.0; n]; n];
        w[0][2] = 1.0;
        w[2][0] = 1.0;
        let mut part = vec![0, 0, 1, 1];
        let valid = |p: &[usize]| p[0] == p[2];
        partition(&w, &mut part, 2, Some(&valid));
        assert_eq!(part[0], part[2]);
    }
}
