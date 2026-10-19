//! Perfect layouts: placing the qubits so that every pair a circuit's
//! two-qubit gates use is coupled, when there is such a placement, so that
//! no SWAP is needed (Qiskit's VF2Layout does the same).
//!
//! It is a subgraph search: the interaction graph (logical qubits, an edge
//! for each pair that meets in a gate) must map into the coupling graph,
//! each logical qubit to a different physical one, each edge to a coupling.
//! Backtracking places the logical qubits in breadth-first order from the
//! busiest, each next to an already placed neighbour, which prunes early; a
//! limit on the steps keeps hopeless searches short.

use super::coupling::CouplingGraph;
use std::collections::BTreeSet;

/// Steps of the search before giving up.
const STEP_LIMIT: usize = 200_000;

/// A layout (`layout[q]`: the physical qubit of logical qubit `q`) under
/// which every pair of `pairs` is coupled, if one is found.
pub fn perfect_layout(
    pairs: &[(usize, usize)],
    logical: usize,
    graph: &CouplingGraph,
) -> Option<Vec<usize>> {
    let mut adjacency: Vec<BTreeSet<usize>> = vec![BTreeSet::new(); logical];
    for &(a, b) in pairs {
        if a != b {
            adjacency[a].insert(b);
            adjacency[b].insert(a);
        }
    }
    let order = placement_order(&adjacency);
    let mut search = Search {
        adjacency: &adjacency,
        graph,
        order: &order,
        layout: vec![usize::MAX; logical],
        used: vec![false; graph.qubits()],
        steps: 0,
    };
    search.place(0).then_some(search.layout)
}

/// Logical qubits breadth-first from the busiest, each component in turn;
/// idle qubits last.
fn placement_order(adjacency: &[BTreeSet<usize>]) -> Vec<usize> {
    let n = adjacency.len();
    let mut seen = vec![false; n];
    let mut order = Vec::with_capacity(n);
    let mut by_degree: Vec<usize> = (0..n).collect();
    by_degree.sort_by_key(|&q| std::cmp::Reverse(adjacency[q].len()));
    for &start in &by_degree {
        if seen[start] || adjacency[start].is_empty() {
            continue;
        }
        seen[start] = true;
        let mut queue = std::collections::VecDeque::from([start]);
        while let Some(q) = queue.pop_front() {
            order.push(q);
            let mut next: Vec<usize> = adjacency[q].iter().copied().filter(|&r| !seen[r]).collect();
            next.sort_by_key(|&r| std::cmp::Reverse(adjacency[r].len()));
            for r in next {
                seen[r] = true;
                queue.push_back(r);
            }
        }
    }
    order.extend((0..n).filter(|&q| !seen[q]));
    order
}

struct Search<'a> {
    adjacency: &'a [BTreeSet<usize>],
    graph: &'a CouplingGraph,
    order: &'a [usize],
    layout: Vec<usize>,
    used: Vec<bool>,
    steps: usize,
}

impl Search<'_> {
    /// Place `order[k..]`; true once every qubit is placed.
    fn place(&mut self, k: usize) -> bool {
        if k == self.order.len() {
            return true;
        }
        self.steps += 1;
        if self.steps > STEP_LIMIT {
            return false;
        }
        let q = self.order[k];
        // Next to an already placed neighbour if it has one, else anywhere.
        let anchor = self.adjacency[q]
            .iter()
            .find(|&&r| self.layout[r] != usize::MAX)
            .map(|&r| self.layout[r]);
        let candidates: Vec<usize> = match anchor {
            Some(p) => self.graph.neighbours(p).to_vec(),
            None => (0..self.graph.qubits()).collect(),
        };
        for p in candidates {
            if self.used[p] || self.graph.neighbours(p).len() < self.adjacency[q].len() {
                continue;
            }
            let fits = self.adjacency[q]
                .iter()
                .all(|&r| self.layout[r] == usize::MAX || self.graph.coupled(p, self.layout[r]));
            if !fits {
                continue;
            }
            self.layout[q] = p;
            self.used[p] = true;
            if self.place(k + 1) {
                return true;
            }
            self.layout[q] = usize::MAX;
            self.used[p] = false;
        }
        false
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_line_fits_a_ring_and_a_triangle_does_not_fit_a_grid() {
        let ring =
            CouplingGraph::new(6, &[(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)]).unwrap();
        let line = [(0, 1), (1, 2), (2, 3), (3, 4)];
        let layout = perfect_layout(&line, 5, &ring).unwrap();
        assert!(
            line.iter()
                .all(|&(a, b)| ring.coupled(layout[a], layout[b]))
        );
        let mut seen = layout.clone();
        seen.sort();
        seen.dedup();
        assert_eq!(seen.len(), 5);
        let grid = CouplingGraph::new(4, &[(0, 1), (1, 3), (3, 2), (2, 0)]).unwrap();
        assert!(perfect_layout(&[(0, 1), (1, 2), (2, 0)], 3, &grid).is_none());
    }
}
