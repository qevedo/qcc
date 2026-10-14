//! Which physical qubits of a chip can interact directly.

use std::collections::VecDeque;

/// The coupling graph of a chip: physical qubits `0..n`, and the pairs that
/// can run a two-qubit gate directly. Shortest-path distances are precomputed.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CouplingGraph {
    n: usize,
    edges: Vec<(usize, usize)>,
    neighbours: Vec<Vec<usize>>,
    /// `distance[a][b]`: edges on a shortest path from `a` to `b`.
    distance: Vec<Vec<u32>>,
}

impl CouplingGraph {
    /// A chip whose qubits are coupled by `edges`, which must connect them all.
    pub fn new(n: usize, edges: &[(usize, usize)]) -> Result<CouplingGraph, String> {
        if n == 0 {
            return Err("a chip needs at least one qubit".into());
        }
        let mut neighbours = vec![Vec::new(); n];
        let mut unique = Vec::new();
        for &(a, b) in edges {
            if a >= n || b >= n || a == b {
                return Err(format!(
                    "the coupling ({a}, {b}) does not join two of the {n} qubits"
                ));
            }
            let edge = (a.min(b), a.max(b));
            if !unique.contains(&edge) {
                unique.push(edge);
                neighbours[a].push(b);
                neighbours[b].push(a);
            }
        }
        unique.sort();
        for list in &mut neighbours {
            list.sort();
        }
        let mut distance = vec![vec![u32::MAX; n]; n];
        for (start, row) in distance.iter_mut().enumerate() {
            row[start] = 0;
            let mut queue = VecDeque::from([start]);
            while let Some(p) = queue.pop_front() {
                for &next in &neighbours[p] {
                    if row[next] == u32::MAX {
                        row[next] = row[p] + 1;
                        queue.push_back(next);
                    }
                }
            }
            if let Some(unreachable) = row.iter().position(|&d| d == u32::MAX) {
                return Err(format!(
                    "the couplings do not connect qubit {unreachable} to qubit {start}"
                ));
            }
        }
        Ok(CouplingGraph {
            n,
            edges: unique,
            neighbours,
            distance,
        })
    }

    /// Every qubit coupled to every other: no routing needed.
    pub fn all_to_all(n: usize) -> CouplingGraph {
        let edges: Vec<_> = (0..n)
            .flat_map(|a| (a + 1..n).map(move |b| (a, b)))
            .collect();
        CouplingGraph::new(n, &edges).expect("connected")
    }

    /// Qubits in a line: 0 – 1 – 2 – …
    pub fn line(n: usize) -> CouplingGraph {
        let edges: Vec<_> = (1..n).map(|p| (p - 1, p)).collect();
        CouplingGraph::new(n, &edges).expect("connected")
    }

    /// Qubits in a ring.
    pub fn ring(n: usize) -> CouplingGraph {
        let mut edges: Vec<_> = (1..n).map(|p| (p - 1, p)).collect();
        if n > 2 {
            edges.push((n - 1, 0));
        }
        CouplingGraph::new(n, &edges).expect("connected")
    }

    /// `rows × columns` qubits in a grid; qubit `r · columns + c` is in row `r`,
    /// column `c`.
    pub fn grid(rows: usize, columns: usize) -> CouplingGraph {
        let mut edges = Vec::new();
        for r in 0..rows {
            for c in 0..columns {
                let p = r * columns + c;
                if c + 1 < columns {
                    edges.push((p, p + 1));
                }
                if r + 1 < rows {
                    edges.push((p, p + columns));
                }
            }
        }
        CouplingGraph::new(rows * columns, &edges).expect("connected")
    }

    pub fn qubits(&self) -> usize {
        self.n
    }

    pub fn edges(&self) -> &[(usize, usize)] {
        &self.edges
    }

    pub fn neighbours(&self, p: usize) -> &[usize] {
        &self.neighbours[p]
    }

    pub fn distance(&self, a: usize, b: usize) -> u32 {
        self.distance[a][b]
    }

    pub fn coupled(&self, a: usize, b: usize) -> bool {
        self.distance[a][b] == 1
    }

    /// A shortest path from `a` to `b`, both included.
    pub fn path(&self, a: usize, b: usize) -> Vec<usize> {
        let mut path = vec![a];
        let mut p = a;
        while p != b {
            p = *self.neighbours[p]
                .iter()
                .find(|&&q| self.distance[q][b] + 1 == self.distance[p][b])
                .expect("connected");
            path.push(p);
        }
        path
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn distances_and_paths() {
        let g = CouplingGraph::grid(3, 3);
        assert_eq!((g.distance(0, 8), g.distance(4, 1)), (4, 1));
        let path = g.path(0, 8);
        assert_eq!((path.len(), path[0], path[4]), (5, 0, 8));
        assert!(path.windows(2).all(|w| g.coupled(w[0], w[1])));
        assert!(
            CouplingGraph::new(3, &[(0, 1)])
                .unwrap_err()
                .contains("do not connect")
        );
        assert_eq!(CouplingGraph::ring(5).distance(0, 3), 2);
    }
}
