//! Communication at fixed places: each link between two cores meets each of
//! their chips at one physical qubit, its *port*.
//!
//! By default ([`physical`](super::physical)), a qubit teleports from any slot
//! of its core to any free slot of another, as if every qubit could reach the
//! link. On hardware, a link ends at a few communication qubits. With ports:
//!
//! - to cross a link, a qubit is first routed (with SWAPs) onto its core's
//!   port of that link; it lands on the other core's port, which is emptied
//!   first by moving the nearest free slot there;
//! - a move between cores that are not linked crosses the links one by one,
//!   from port to port through the cores in between;
//! - the mapper keeps one slot free in every core, so a core can always take
//!   a qubit in before another leaves, and the moves before a timeslice never
//!   wait for each other in a cycle (as TeleSABRE does, Russo et al., QCE 2025).
//!
//! The moves before a timeslice run one after another, each when its next
//! core has room, crossing links on any shortest path. When every next core
//! is full (two qubits crossing between full cores, for instance), two
//! qubits exchange places across a link: a teleportation each way at once.
//! If no pair is moving that way, one of them moves back later, and those
//! teleportations are counted on top of the mapping's.

use super::architecture::Architecture;
use super::physical::Builder;
use crate::ir::Instruction;
use crate::routing::CouplingGraph;

/// The port of every link in each of its two cores.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Ports {
    /// `slots[l]`: for link `l` of [`Architecture::links`], `(a, b)` with the
    /// smaller core first, the slot of its port in core `a` and in core `b`.
    slots: Vec<(usize, usize)>,
}

impl Ports {
    /// Ports at the given slots, one pair per link, in the order of
    /// [`Architecture::links`].
    pub fn new(slots: Vec<(usize, usize)>) -> Ports {
        Ports { slots }
    }

    /// Ports at the most central free slots of each core's chip (those with
    /// the smallest total distance to every slot), one per link, so qubits
    /// have least far to go.
    pub fn central(arch: &Architecture, chip: &CouplingGraph) -> Ports {
        let k = chip.qubits();
        let mut by_centrality: Vec<usize> = (0..k).collect();
        by_centrality.sort_by_key(|&p| ((0..k).map(|o| chip.distance(p, o)).sum::<u32>(), p));
        let mut next = vec![0; arch.cores];
        let mut take = |core: usize| {
            let slot = by_centrality[next[core] % k];
            next[core] += 1;
            slot
        };
        let slots = arch
            .links()
            .iter()
            .map(|&(a, b)| (take(a), take(b)))
            .collect();
        Ports { slots }
    }

    /// Check that there is one port pair per link, on the chip, with no slot
    /// serving two links of the same core.
    pub fn check(&self, arch: &Architecture, chip: &CouplingGraph) -> Result<(), String> {
        if self.slots.len() != arch.links().len() {
            return Err(format!(
                "{} port pairs for {} links",
                self.slots.len(),
                arch.links().len()
            ));
        }
        let mut used = vec![Vec::new(); arch.cores];
        for (&(a, b), &(sa, sb)) in arch.links().iter().zip(&self.slots) {
            for (core, slot) in [(a, sa), (b, sb)] {
                if slot >= chip.qubits() {
                    return Err(format!(
                        "port slot {slot} of core {core} is beyond its {} qubits",
                        chip.qubits()
                    ));
                }
                if used[core].contains(&slot) {
                    return Err(format!(
                        "slot {slot} of core {core} is the port of two links"
                    ));
                }
                used[core].push(slot);
            }
        }
        Ok(())
    }

    /// The mean number of SWAPs from a slot of a core to one of its ports,
    /// over every slot and port.
    pub fn mean_distance(&self, chip: &CouplingGraph) -> f64 {
        let k = chip.qubits();
        let total: u32 = self
            .slots
            .iter()
            .flat_map(|&(a, b)| [a, b])
            .map(|port| (0..k).map(|s| chip.distance(s, port)).sum::<u32>())
            .sum();
        f64::from(total) / (2 * self.slots.len() * k).max(1) as f64
    }

    /// The slot of `core`'s port on its link to `other`.
    pub fn slot(&self, arch: &Architecture, core: usize, other: usize) -> usize {
        let link = (core.min(other), core.max(other));
        let l = arch
            .links()
            .binary_search(&link)
            .expect("the cores are linked");
        if core < other {
            self.slots[l].0
        } else {
            self.slots[l].1
        }
    }
}

impl Builder<'_> {
    /// Before the first timeslice: put each core's qubits that stay where
    /// SABRE's layout search wants them for the first epoch, and the ones
    /// that leave on the free slots closest to the port they leave by.
    pub(super) fn place_initial(
        &mut self,
        work: &[&Instruction],
        assignment: &[usize],
        initial: &[usize],
        ports: &Ports,
    ) -> Result<(), String> {
        let k = self.k;
        for core in 0..self.arch.cores {
            let stay: Vec<usize> = Self::members(initial, core)
                .into_iter()
                .filter(|&q| assignment[q] == core)
                .collect();
            let local = super::physical::local_circuit(work, &stay);
            let layout = if local.is_empty() {
                (0..stay.len()).collect()
            } else {
                crate::routing::layout_and_route(
                    &local,
                    stay.len(),
                    self.chip,
                    &self.options.sabre,
                )?
                .initial_layout
            };
            let mut free = vec![true; k];
            for (&q, &slot) in stay.iter().zip(&layout) {
                self.slot_of[q] = core * k + slot;
                free[slot] = false;
            }
            for q in Self::members(initial, core) {
                if assignment[q] == core {
                    continue;
                }
                let next = self.next_core(core, assignment[q]);
                let port = ports.slot(self.arch, core, next);
                let slot = (0..k)
                    .filter(|&s| free[s])
                    .min_by_key(|&s| (self.chip.distance(s, port), s))
                    .ok_or("the initial assignment overfills a core")?;
                free[slot] = false;
                self.slot_of[q] = core * k + slot;
            }
        }
        for (q, &p) in self.slot_of.iter().enumerate() {
            self.occupant[p] = Some(q);
        }
        self.initial_layout = self.slot_of.clone();
        Ok(())
    }

    /// Move the qubits whose core changes in `assignment` through the ports,
    /// one link at a time, each when the core it enters has a free slot.
    pub(super) fn teleport_through_ports(
        &mut self,
        assignment: &[usize],
        ports: &Ports,
    ) -> Result<(), String> {
        let k = self.k;
        let mut exchanges = 0;
        loop {
            let pending: Vec<usize> = (0..assignment.len())
                .filter(|&q| self.slot_of[q] / k != assignment[q])
                .collect();
            if pending.is_empty() {
                return Ok(());
            }
            // The first move with a next core, on a shortest path, that has room.
            let ready = pending.iter().find_map(|&q| {
                let here = self.slot_of[q] / k;
                self.next_cores(here, assignment[q])
                    .find(|&next| self.free_slots(next) > 0)
                    .map(|next| (q, here, next))
            });
            match ready {
                Some((q, here, next)) => self.hop(q, here, next, ports),
                None => {
                    // Every next core is full: two qubits exchange places
                    // across a link instead.
                    exchanges += 1;
                    if exchanges > 10 * assignment.len() {
                        return Err("qubits keep exchanging places between full cores".into());
                    }
                    let (q, here, next, r) = self.exchange_pair(&pending, assignment);
                    self.exchange(q, here, r, next, ports);
                }
            }
        }
    }

    /// When every move's next core is full: a moving qubit and a qubit of its
    /// next core to exchange with. Preferably one moving the other way across
    /// the same link (then neither moves more than the mapping says); else
    /// one that moves anywhere, else any (which then moves back, costing
    /// teleportations the mapping did not count).
    fn exchange_pair(
        &self,
        pending: &[usize],
        assignment: &[usize],
    ) -> (usize, usize, usize, usize) {
        let k = self.k;
        let core = |q: usize| self.slot_of[q] / k;
        for &q in pending {
            for next in self.next_cores(core(q), assignment[q]) {
                if let Some(&r) = pending.iter().find(|&&r| {
                    core(r) == next && self.next_cores(next, assignment[r]).any(|c| c == core(q))
                }) {
                    return (q, core(q), next, r);
                }
            }
        }
        let q = pending[0];
        let next = self.next_core(core(q), assignment[q]);
        let r = pending
            .iter()
            .copied()
            .find(|&r| core(r) == next)
            .or_else(|| (next * k..(next + 1) * k).find_map(|p| self.occupant[p]))
            .expect("a full core has qubits");
        (q, core(q), next, r)
    }

    /// Exchange `q` (in core `a`) and `r` (in the linked core `b`) across their
    /// link: both onto their ports, then a teleportation each way at once, as
    /// with a communication qubit at each end. Barriers on the first port keep
    /// the pair apart from other teleportations, which would otherwise read as
    /// happening at the same time.
    fn exchange(&mut self, q: usize, a: usize, r: usize, b: usize, ports: &Ports) {
        let k = self.k;
        let pa = a * k + ports.slot(self.arch, a, b);
        let pb = b * k + ports.slot(self.arch, b, a);
        self.walk(q, pa);
        self.walk(r, pb);
        let barrier = Instruction::gate("barrier", vec![pa], vec![]);
        self.out.push(barrier.clone());
        self.out
            .push(Instruction::gate("teleport", vec![pa, pb], vec![]));
        self.out
            .push(Instruction::gate("teleport", vec![pb, pa], vec![]));
        self.out.push(barrier);
        self.occupant[pa] = Some(r);
        self.occupant[pb] = Some(q);
        self.slot_of[q] = pb;
        self.slot_of[r] = pa;
        self.teleportations += 2;
    }

    /// Move `q` onto physical qubit `p` of its core, with SWAPs along a
    /// shortest path.
    fn walk(&mut self, q: usize, p: usize) {
        let base = p / self.k * self.k;
        let path = self.chip.path(self.slot_of[q] - base, p - base);
        for w in path.windows(2) {
            self.swap(base + w[0], base + w[1]);
        }
    }

    /// Teleport `q` from core `a` to the linked core `b`: onto `a`'s port,
    /// across, onto `b`'s port (emptied first).
    fn hop(&mut self, q: usize, a: usize, b: usize, ports: &Ports) {
        let k = self.k;
        let from = a * k + ports.slot(self.arch, a, b);
        let to = b * k + ports.slot(self.arch, b, a);
        self.walk(q, from);
        // Bring the nearest free slot of b onto its port.
        if self.occupant[to].is_some() {
            let port = to - b * k;
            let hole = (0..k)
                .filter(|&s| self.occupant[b * k + s].is_none())
                .min_by_key(|&s| (self.chip.distance(s, port), s))
                .expect("the core has a free slot");
            let path = self.chip.path(hole, port);
            for w in path.windows(2) {
                self.swap(b * k + w[0], b * k + w[1]);
            }
        }
        self.out
            .push(Instruction::gate("teleport", vec![from, to], vec![]));
        self.occupant[from] = None;
        self.occupant[to] = Some(q);
        self.slot_of[q] = to;
        self.teleportations += 1;
    }

    /// A SWAP inside a core (one of the qubits may be empty).
    fn swap(&mut self, x: usize, y: usize) {
        self.out.push(Instruction::gate("swap", vec![x, y], vec![]));
        let (qx, qy) = (self.occupant[x], self.occupant[y]);
        self.occupant[x] = qy;
        self.occupant[y] = qx;
        if let Some(q) = qx {
            self.slot_of[q] = y;
        }
        if let Some(q) = qy {
            self.slot_of[q] = x;
        }
        self.swaps += 1;
    }

    /// The cores next to `a` on a shortest path of links to `b`.
    fn next_cores(&self, a: usize, b: usize) -> impl Iterator<Item = usize> + '_ {
        (0..self.arch.cores).filter(move |&n| {
            self.arch.hops(a, n) == 1 && self.arch.hops(n, b) + 1 == self.arch.hops(a, b)
        })
    }

    /// The first of [`Self::next_cores`].
    fn next_core(&self, a: usize, b: usize) -> usize {
        self.next_cores(a, b).next().unwrap_or(b)
    }

    fn free_slots(&self, core: usize) -> usize {
        let k = self.k;
        (core * k..(core + 1) * k)
            .filter(|&p| self.occupant[p].is_none())
            .count()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::multicore::{Method, PhysicalOptions, Timing, map_and_route};
    use crate::rng::Rng;

    /// Every method, on many random machines with ports and circuits that
    /// fill them, moves every qubit through the ports.
    #[test]
    fn moves_never_get_stuck() {
        let mut rng = Rng::new(11);
        let (mut extra, mut total) = (0, 0);
        for trial in 0..60 {
            let cores = 2 + rng.below(5);
            let k = 3 + rng.below(4);
            let arch = match trial % 3 {
                0 => Architecture::line(cores, k),
                1 => Architecture::ring(cores, k),
                _ => Architecture::grid(2, cores.div_ceil(2), k),
            }
            .unwrap();
            let chip = CouplingGraph::line(k);
            let qubits = arch.cores * (k - 1) - rng.below(2);
            let circuit: Vec<Instruction> = (0..8 * qubits)
                .map(|_| {
                    let a = rng.below(qubits);
                    let b = (a + 1 + rng.below(qubits - 1)) % qubits;
                    Instruction::gate("cx", vec![a, b], vec![])
                })
                .collect();
            let options = PhysicalOptions {
                ports: Some(Ports::central(&arch, &chip)),
                sabre: crate::routing::SabreOptions {
                    trials: 1,
                    ..Default::default()
                },
                ..PhysicalOptions::default()
            };
            for method in [Method::Hqa, Method::FgpRoee, Method::Naive] {
                let (mapping, physical) = map_and_route(
                    &circuit,
                    qubits,
                    &arch,
                    &chip,
                    method,
                    &Timing::default(),
                    &options,
                )
                .unwrap_or_else(|e| panic!("{method:?} on {trial}: {e}"));
                // Exchanges may add teleportations, but rarely.
                assert!(physical.teleportations >= mapping.communications());
                extra += physical.teleportations - mapping.communications();
                total += mapping.communications();
            }
        }
        assert!(
            extra * 100 <= total,
            "{extra} extra teleportations over {total}"
        );
    }

    #[test]
    fn central_ports_are_distinct_and_checked() {
        let arch = Architecture::grid(2, 2, 9).unwrap();
        let chip = CouplingGraph::grid(3, 3);
        let ports = Ports::central(&arch, &chip);
        ports.check(&arch, &chip).unwrap();
        // Each core has two links: the centre and then an edge midpoint.
        assert_eq!(ports.slot(&arch, 0, 1), 4);
        assert_ne!(ports.slot(&arch, 0, 2), 4);
        assert!(Ports::new(vec![(0, 0); 3]).check(&arch, &chip).is_err());
        assert!(
            Ports::new(vec![(0, 0), (0, 1), (1, 1), (2, 2)])
                .check(&arch, &chip)
                .unwrap_err()
                .contains("two links")
        );
    }
}
