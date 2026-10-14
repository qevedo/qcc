//! The hardware model: a modular quantum computer made of several cores.

/// A multi-core quantum computer: `cores` cores of `qubits_per_core` physical
/// qubits each, connected by quantum links.
///
/// Every core is all-to-all inside (any two of its qubits can interact
/// directly). A two-qubit gate can only run when both of its qubits are in the
/// same core, and moving a qubit's state over one link is one teleportation,
/// one *non-local communication*; moving it between cores that are not linked
/// takes one per link on the shortest path (`hops`).
///
/// Escofet et al. link all cores to all ([`Architecture::new`]); lines, rings,
/// grids and any connected set of links are also available.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Architecture {
    pub cores: usize,
    pub qubits_per_core: usize,
    /// The links between cores, each pair once with the smaller core first.
    links: Vec<(usize, usize)>,
    /// `hops[a][b]`: links on a shortest path from core `a` to core `b`.
    hops: Vec<Vec<usize>>,
}

impl Architecture {
    /// All cores linked to all, as in the papers.
    pub fn new(cores: usize, qubits_per_core: usize) -> Result<Architecture, String> {
        let links = (0..cores)
            .flat_map(|a| (a + 1..cores).map(move |b| (a, b)))
            .collect::<Vec<_>>();
        Architecture::with_links(cores, qubits_per_core, &links)
    }

    /// Cores linked by `links` (pairs of core indices), which must connect them all.
    pub fn with_links(
        cores: usize,
        qubits_per_core: usize,
        links: &[(usize, usize)],
    ) -> Result<Architecture, String> {
        if cores < 1 {
            return Err("an architecture needs at least one core".into());
        }
        if qubits_per_core < 2 {
            return Err("each core needs at least two qubits to run a two-qubit gate".into());
        }
        let mut neighbours = vec![Vec::new(); cores];
        let mut unique = Vec::new();
        for &(a, b) in links {
            if a >= cores || b >= cores || a == b {
                return Err(format!(
                    "the link ({a}, {b}) does not join two of the {cores} cores"
                ));
            }
            let link = (a.min(b), a.max(b));
            if !unique.contains(&link) {
                unique.push(link);
                neighbours[a].push(b);
                neighbours[b].push(a);
            }
        }
        unique.sort();
        // Shortest paths by breadth-first search from every core.
        let mut hops = vec![vec![usize::MAX; cores]; cores];
        for (start, row) in hops.iter_mut().enumerate() {
            row[start] = 0;
            let mut queue = std::collections::VecDeque::from([start]);
            while let Some(core) = queue.pop_front() {
                for &next in &neighbours[core] {
                    if row[next] == usize::MAX {
                        row[next] = row[core] + 1;
                        queue.push_back(next);
                    }
                }
            }
            if let Some(unreachable) = row.iter().position(|&h| h == usize::MAX) {
                return Err(format!(
                    "the links do not connect core {unreachable} to core {start}"
                ));
            }
        }
        Ok(Architecture {
            cores,
            qubits_per_core,
            links: unique,
            hops,
        })
    }

    /// Cores in a line: 0 – 1 – 2 – …
    pub fn line(cores: usize, qubits_per_core: usize) -> Result<Architecture, String> {
        let links: Vec<_> = (1..cores).map(|c| (c - 1, c)).collect();
        Architecture::with_links(cores, qubits_per_core, &links)
    }

    /// Cores in a ring: a line whose ends are linked.
    pub fn ring(cores: usize, qubits_per_core: usize) -> Result<Architecture, String> {
        let mut links: Vec<_> = (1..cores).map(|c| (c - 1, c)).collect();
        if cores > 2 {
            links.push((cores - 1, 0));
        }
        Architecture::with_links(cores, qubits_per_core, &links)
    }

    /// `rows × columns` cores in a grid, each linked to its four neighbours.
    /// Core `r · columns + c` is in row `r`, column `c`.
    pub fn grid(
        rows: usize,
        columns: usize,
        qubits_per_core: usize,
    ) -> Result<Architecture, String> {
        let mut links = Vec::new();
        for r in 0..rows {
            for c in 0..columns {
                let core = r * columns + c;
                if c + 1 < columns {
                    links.push((core, core + 1));
                }
                if r + 1 < rows {
                    links.push((core, core + columns));
                }
            }
        }
        Architecture::with_links(rows * columns, qubits_per_core, &links)
    }

    /// The total number of physical qubits.
    pub fn capacity(&self) -> usize {
        self.cores * self.qubits_per_core
    }

    /// Teleportations needed to move a qubit from core `a` to core `b`.
    pub fn hops(&self, a: usize, b: usize) -> usize {
        self.hops[a][b]
    }

    /// The links between cores, each once with the smaller core first.
    pub fn links(&self) -> &[(usize, usize)] {
        &self.links
    }

    /// Whether every core is linked to every other.
    pub fn is_all_to_all(&self) -> bool {
        self.links.len() == self.cores * (self.cores - 1) / 2
    }

    /// Check that `qubits` circuit qubits fit.
    pub fn check_fits(&self, qubits: usize) -> Result<(), String> {
        if qubits > self.capacity() {
            return Err(format!(
                "{qubits} qubits do not fit in {} cores of {} qubits",
                self.cores, self.qubits_per_core
            ));
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hops_follow_the_links() {
        let all = Architecture::new(4, 2).unwrap();
        assert!(all.is_all_to_all());
        assert_eq!(all.hops(0, 3), 1);
        let line = Architecture::line(4, 2).unwrap();
        assert_eq!((line.hops(0, 3), line.hops(2, 1)), (3, 1));
        let ring = Architecture::ring(6, 2).unwrap();
        assert_eq!((ring.hops(0, 5), ring.hops(0, 3)), (1, 3));
        let grid = Architecture::grid(2, 3, 2).unwrap();
        assert_eq!((grid.hops(0, 5), grid.hops(1, 4)), (3, 1));
        assert!(!grid.is_all_to_all());
    }

    #[test]
    fn links_must_connect_every_core() {
        assert!(
            Architecture::with_links(3, 2, &[(0, 1)])
                .unwrap_err()
                .contains("do not connect")
        );
        assert!(
            Architecture::with_links(3, 2, &[(0, 3)])
                .unwrap_err()
                .contains("does not join")
        );
        let a = Architecture::with_links(3, 2, &[(1, 0), (0, 1), (2, 1)]).unwrap();
        assert_eq!(a.links(), &[(0, 1), (1, 2)]);
    }
}
