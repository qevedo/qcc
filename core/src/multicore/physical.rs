//! From a mapping to a circuit that runs: routing inside each core and
//! teleportations between cores.
//!
//! A [`Mapping`] says which core each qubit is in during each timeslice, as if
//! every core were all-to-all inside. Real cores are not: their qubits are
//! coupled in a line, a grid or another graph, the *chip* of each core. This
//! module turns a mapping into instructions on physical qubits:
//!
//! 1. The timeslices are grouped into **epochs**, runs of slices in which no
//!    qubit changes core.
//! 2. Inside each epoch, each core runs its part of the circuit on its own
//!    chip, routed with SABRE ([`crate::routing`]), which inserts SWAPs. Each
//!    core starts where the previous epoch left its qubits, and SABRE also
//!    looks at the gates its qubits have next, in the following epochs.
//! 3. Between epochs, the qubits that change core are teleported. A qubit
//!    arriving in a core takes a free slot close to the qubit it will interact
//!    with next.
//!
//! The first epoch's layout inside each core is found with SABRE's layout
//! search; the initial placement is free, as in the mappers.
//!
//! Physical qubit `core · qubits_per_core + slot` is slot `slot` of core
//! `core`. Teleportations are `"teleport"` instructions on (source,
//! destination): the qubit's state moves from the source, which is left empty,
//! to the destination, which was empty or is being left by another
//! teleportation before the same timeslice. The teleportations before a
//! timeslice happen at once, as with dedicated communication qubits in each
//! core; they are listed so that each destination is already free whenever
//! possible, which only fails when qubits swap places between full cores.
//! A teleportation between cores that are not linked is one instruction, but
//! counts one teleportation per link it crosses.
//!
//! With [`Ports`] ([`PhysicalOptions::ports`]), communication happens at fixed
//! places instead, as on hardware: see [`super::ports`].

use super::architecture::Architecture;
use super::assignment::Mapping;
use super::ports::Ports;
use super::timeslices::Timeslices;
use crate::ir::Instruction;
use crate::routing::{CouplingGraph, SabreOptions, layout_and_route, route, route_ahead};

/// How to route inside the cores.
#[derive(Clone, Debug)]
pub struct PhysicalOptions {
    pub sabre: SabreOptions,
    /// Put a qubit arriving in a core next to its next partner (or where its
    /// partner can land next to it). Without it, it takes the first free slot.
    pub partner_slots: bool,
    /// How many timeslices past an epoch SABRE looks at (0: none). The default,
    /// 4, was chosen on the benchmarks; 1 to 64 all give totals within 2%.
    pub look_ahead: usize,
    /// Where each link meets each core's chip ([`Ports`]). With ports, a qubit
    /// must be routed onto its core's port of a link to teleport, and lands on
    /// the other core's port; the mapper then keeps a slot free in every core.
    /// Without (the default), a qubit teleports from any slot to any free one.
    pub ports: Option<Ports>,
}

impl Default for PhysicalOptions {
    fn default() -> PhysicalOptions {
        PhysicalOptions {
            sabre: SabreOptions::default(),
            partner_slots: true,
            look_ahead: 4,
            ports: None,
        }
    }
}

/// A multi-core circuit on physical qubits, ready to run.
#[derive(Clone, Debug)]
pub struct PhysicalCircuit {
    /// The instructions on physical qubits (`core · qubits_per_core + slot`),
    /// with SWAPs (`"swap"`) inside cores and teleportations (`"teleport"`)
    /// between them. Barriers are dropped.
    pub instructions: Vec<Instruction>,
    /// `initial_layout[q]`: the physical qubit of circuit qubit `q` at the start.
    pub initial_layout: Vec<usize>,
    /// `final_layout[q]`: where circuit qubit `q` ends up.
    pub final_layout: Vec<usize>,
    /// SWAPs inserted inside cores.
    pub swaps: usize,
    /// Teleportations, one per link crossed: the mapping's communications.
    pub teleportations: usize,
}

/// Turn `circuit` (on qubits `0..slices.qubits()`), mapped onto cores by
/// `mapping` for its `slices`, into a circuit on physical qubits, where the
/// qubits of each core are coupled as in `chip`.
///
/// `slices` must be the timeslices the mapping was made for, built from the
/// same circuit (as [`map_and_route`] does).
pub fn physical_circuit(
    circuit: &[Instruction],
    mapping: &Mapping,
    slices: &Timeslices,
    chip: &CouplingGraph,
    options: &PhysicalOptions,
) -> Result<PhysicalCircuit, String> {
    let arch = &mapping.architecture;
    // With ports, the mapper leaves a slot free in every core.
    let reserved = usize::from(options.ports.is_some());
    if chip.qubits() < arch.qubits_per_core + reserved {
        return Err(format!(
            "the chip has {} qubits, the mapping fills {} per core{}",
            chip.qubits(),
            arch.qubits_per_core,
            if reserved == 1 {
                " and ports need one free"
            } else {
                ""
            }
        ));
    }
    if let Some(ports) = &options.ports {
        ports.check(arch, chip)?;
    }
    if mapping.slices.len() != slices.len() {
        return Err(format!(
            "{} assignments for {} timeslices",
            mapping.slices.len(),
            slices.len()
        ));
    }
    let stages = stages(circuit, slices)?;
    let epochs = epochs(mapping);
    let mut builder = Builder::new(arch, chip, options, slices.qubits());

    // Instructions per epoch, in program order.
    let mut epoch_of_stage = vec![0; slices.len().max(1)];
    for (e, &(start, end)) in epochs.iter().enumerate() {
        epoch_of_stage[start..end.max(start + 1)].fill(e);
    }
    let mut work: Vec<Vec<&Instruction>> = vec![Vec::new(); epochs.len()];
    for (inst, stage) in circuit.iter().zip(&stages) {
        if let Some(t) = stage {
            work[epoch_of_stage[*t]].push(inst);
        }
    }

    let assignment_at = |t: usize| mapping.slices.get(t).unwrap_or(&mapping.initial);
    for (e, &(start, end)) in epochs.iter().enumerate() {
        let assignment = assignment_at(start);
        match (e, &options.ports) {
            (0, None) => builder.first_epoch(&work[0], assignment, &mapping.initial)?,
            (0, Some(ports)) => {
                builder.place_initial(&work[0], assignment, &mapping.initial, ports)?;
                builder.teleport_through_ports(assignment, ports)?;
                builder.epoch(&work[0], assignment, mapping, slices, end)?;
            }
            (_, None) => {
                builder.teleport(assignment, slices, start);
                builder.epoch(&work[e], assignment, mapping, slices, end)?;
            }
            (_, Some(ports)) => {
                builder.teleport_through_ports(assignment, ports)?;
                builder.epoch(&work[e], assignment, mapping, slices, end)?;
            }
        }
    }
    Ok(PhysicalCircuit {
        instructions: builder.out,
        initial_layout: builder.initial_layout,
        final_layout: builder.slot_of,
        swaps: builder.swaps,
        teleportations: builder.teleportations,
    })
}

/// Map `circuit` onto `arch` with `method` and route it inside the cores,
/// whose qubits are coupled as in `chip`. Gates on three or more qubits must
/// be decomposed first. With [`Method::HqaLatency`](super::Method), every
/// latency weight's mapping is routed and the one whose circuit runs fastest
/// under `timing` is kept; through ports, the weights are also tried with a
/// teleportation made dearer by the mean SWAPs to and from the ports.
pub fn map_and_route(
    circuit: &[Instruction],
    qubits: usize,
    arch: &Architecture,
    chip: &CouplingGraph,
    method: super::Method,
    timing: &super::Timing,
    options: &PhysicalOptions,
) -> Result<(Mapping, PhysicalCircuit), String> {
    if chip.qubits() != arch.qubits_per_core {
        return Err(format!(
            "the chip has {} qubits, the cores {}",
            chip.qubits(),
            arch.qubits_per_core
        ));
    }
    // With ports, the mapper plans for one slot fewer per core.
    let planning = match options.ports {
        Some(_) => Architecture::with_links(arch.cores, arch.qubits_per_core - 1, arch.links())?,
        None => arch.clone(),
    };
    let slices = Timeslices::from_instructions(qubits, circuit)?.fitted_to(&planning);
    if method == super::Method::HqaLatency {
        // Route every candidate and keep the one that runs fastest, SWAPs
        // included (the mapper's own estimate leaves them out).
        let initial = super::initial_placement(&slices, &planning);
        let mut candidates = super::latency_candidates(&slices, &planning, &initial, timing)?;
        // Through ports, a teleportation also costs the SWAPs to its port and
        // from the other one: also plan for that, on average.
        if let Some(ports) = &options.ports {
            let planned = super::Timing {
                teleport: timing.teleport + 2.0 * ports.mean_distance(chip) * timing.swap,
                ..*timing
            };
            candidates.extend(super::latency_candidates(
                &slices, &planning, &initial, &planned,
            )?);
        }
        let mut routed = Vec::new();
        for mapping in candidates {
            let physical = physical_circuit(circuit, &mapping, &slices, chip, options)?;
            let duration = super::schedule(&physical.instructions, arch, timing)?.duration;
            routed.push((duration, mapping, physical));
        }
        let (_, mapping, physical) = super::fastest(routed, |(d, m, _)| (*d, m.communications()));
        return Ok((mapping, physical));
    }
    let mapping = super::map_timed(&slices, &planning, method, timing)?;
    let physical = physical_circuit(circuit, &mapping, &slices, chip, options)?;
    Ok((mapping, physical))
}

/// The timeslice each instruction belongs to: its own slice for a two-qubit
/// gate; for anything on one qubit, the slice of the qubit's last two-qubit
/// gate before it (0 if none). Barriers get `None` and are dropped.
fn stages(circuit: &[Instruction], slices: &Timeslices) -> Result<Vec<Option<usize>>, String> {
    let n = slices.qubits();
    // `next[q]`: the first slice where `q`'s next two-qubit gate can be.
    let mut next = vec![0usize; n];
    let mut out = Vec::with_capacity(circuit.len());
    for inst in circuit {
        if let Some(&q) = inst.qubits.iter().find(|&&q| q >= n) {
            return Err(format!(
                "'{}' uses qubit {q}, beyond the {n} of the circuit",
                inst.name
            ));
        }
        let stage = match inst.qubits.as_slice() {
            _ if inst.name == "barrier" => None,
            [] => None,
            &[q] => Some(next[q].saturating_sub(1)),
            &[a, b] if inst.is_gate() => {
                let t = (next[a].max(next[b])..slices.len())
                    .find(|&t| slices.partner(t, a) == Some(b))
                    .ok_or_else(|| {
                        format!(
                            "the gate '{}' on ({a}, {b}) is not in the timeslices",
                            inst.name
                        )
                    })?;
                next[a] = t + 1;
                next[b] = t + 1;
                Some(t)
            }
            qubits => {
                return Err(format!(
                    "'{}' acts on {} qubits: decompose the circuit into one- and two-qubit gates first",
                    inst.name,
                    qubits.len()
                ));
            }
        };
        out.push(stage);
    }
    Ok(out)
}

/// Runs of consecutive timeslices with the same assignment, as `start..end`.
/// A circuit with no two-qubit gates is one empty epoch.
fn epochs(mapping: &Mapping) -> Vec<(usize, usize)> {
    let mut out: Vec<(usize, usize)> = Vec::new();
    for t in 0..mapping.slices.len() {
        match out.last_mut() {
            Some(last) if mapping.slices[t] == mapping.slices[t - 1] => last.1 = t + 1,
            _ => out.push((t, t + 1)),
        }
    }
    if out.is_empty() {
        out.push((0, 0));
    }
    out
}

/// The physical circuit as it is built, epoch by epoch.
pub(super) struct Builder<'a> {
    pub(super) arch: &'a Architecture,
    pub(super) chip: &'a CouplingGraph,
    pub(super) options: &'a PhysicalOptions,
    /// Physical qubits per core: the chip's.
    pub(super) k: usize,
    /// `slot_of[q]`: the physical qubit of circuit qubit `q` now.
    pub(super) slot_of: Vec<usize>,
    /// `occupant[p]`: the circuit qubit on physical qubit `p`, if any.
    pub(super) occupant: Vec<Option<usize>>,
    pub(super) initial_layout: Vec<usize>,
    pub(super) out: Vec<Instruction>,
    pub(super) swaps: usize,
    pub(super) teleportations: usize,
}

impl<'a> Builder<'a> {
    fn new(
        arch: &'a Architecture,
        chip: &'a CouplingGraph,
        options: &'a PhysicalOptions,
        qubits: usize,
    ) -> Builder<'a> {
        Builder {
            arch,
            chip,
            options,
            k: chip.qubits(),
            slot_of: vec![0; qubits],
            occupant: vec![None; arch.cores * chip.qubits()],
            initial_layout: Vec::new(),
            out: Vec::new(),
            swaps: 0,
            teleportations: 0,
        }
    }

    /// The circuit qubits in `core` under `assignment`.
    pub(super) fn members(assignment: &[usize], core: usize) -> Vec<usize> {
        (0..assignment.len())
            .filter(|&q| assignment[q] == core)
            .collect()
    }

    /// The first epoch: find each core's layout with SABRE's layout search,
    /// then work out where the qubits that move before the first timeslice
    /// start, in their initial cores.
    fn first_epoch(
        &mut self,
        work: &[&Instruction],
        assignment: &[usize],
        initial: &[usize],
    ) -> Result<(), String> {
        let k = self.k;
        let mut routed_cores = Vec::new();
        for core in 0..self.arch.cores {
            let members = Self::members(assignment, core);
            let local = local_circuit(work, &members);
            let routed = if local.is_empty() {
                let in_order: Vec<usize> = (0..members.len()).collect();
                route(
                    &local,
                    members.len(),
                    self.chip,
                    &in_order,
                    &self.options.sabre,
                )?
            } else {
                layout_and_route(&local, members.len(), self.chip, &self.options.sabre)?
            };
            for (i, &q) in members.iter().enumerate() {
                self.slot_of[q] = core * k + routed.initial_layout[i];
            }
            routed_cores.push((core, members, routed));
        }
        // Before the first timeslice, the qubits that stay keep their slots;
        // the ones that move start in a slot their initial core does not use.
        let mut start = self.slot_of.clone();
        let mut used = vec![false; self.occupant.len()];
        for q in 0..start.len() {
            if initial[q] == assignment[q] {
                used[start[q]] = true;
            }
        }
        for q in 0..start.len() {
            if initial[q] != assignment[q] {
                let p = (initial[q] * k..(initial[q] + 1) * k)
                    .find(|&p| !used[p])
                    .expect("the initial assignment fits the cores");
                used[p] = true;
                start[q] = p;
            }
        }
        self.initial_layout = start.clone();
        let moves: Vec<(usize, usize, usize)> = (0..start.len())
            .filter(|&q| start[q] != self.slot_of[q])
            .map(|q| (q, start[q], self.slot_of[q]))
            .collect();
        self.slot_of = start;
        for (q, &p) in self.slot_of.iter().enumerate() {
            self.occupant[p] = Some(q);
        }
        self.apply_moves(moves);
        for (core, members, routed) in routed_cores {
            self.emit_routed(core, &members, routed);
        }
        Ok(())
    }

    /// A later epoch: route each core from where its qubits are.
    pub(super) fn epoch(
        &mut self,
        work: &[&Instruction],
        assignment: &[usize],
        mapping: &Mapping,
        slices: &Timeslices,
        end: usize,
    ) -> Result<(), String> {
        let k = self.k;
        for core in 0..self.arch.cores {
            let members = Self::members(assignment, core);
            let local = local_circuit(work, &members);
            if local.is_empty() {
                continue;
            }
            let layout: Vec<usize> = members
                .iter()
                .map(|&q| self.slot_of[q] - core * k)
                .collect();
            let ahead = self.ahead(&members, core, mapping, slices, end);
            let routed = route_ahead(
                &local,
                members.len(),
                self.chip,
                &layout,
                &ahead,
                &self.options.sabre,
            )?;
            self.emit_routed(core, &members, routed);
        }
        Ok(())
    }

    /// The next gates, from timeslice `end` on, between `members` of `core`
    /// that are still in it then, relabelled to `0..members.len()`: SABRE
    /// looks at them to leave the core ready for the next epochs.
    fn ahead(
        &self,
        members: &[usize],
        core: usize,
        mapping: &Mapping,
        slices: &Timeslices,
        end: usize,
    ) -> Vec<(usize, usize)> {
        let mut local = vec![None; mapping.initial.len()];
        for (i, &q) in members.iter().enumerate() {
            local[q] = Some(i);
        }
        let mut out = Vec::new();
        for t in end..slices.len().min(end + self.options.look_ahead) {
            for &(a, b) in slices.slice(t) {
                if let (Some(la), Some(lb)) = (local[a], local[b])
                    && mapping.slices[t][a] == core
                    && mapping.slices[t][b] == core
                {
                    out.push((la, lb));
                    if out.len() == self.options.sabre.extended_set {
                        return out;
                    }
                }
            }
        }
        out
    }

    /// Emit a core's routed circuit on physical qubits and follow its SWAPs.
    fn emit_routed(&mut self, core: usize, members: &[usize], routed: crate::routing::Routed) {
        let base = core * self.k;
        for inst in routed.instructions {
            self.out.push(Instruction {
                qubits: inst.qubits.iter().map(|&p| base + p).collect(),
                ..inst
            });
        }
        self.swaps += routed.swaps;
        for p in base..base + self.k {
            self.occupant[p] = None;
        }
        for (i, &q) in members.iter().enumerate() {
            let p = base + routed.final_layout[i];
            self.slot_of[q] = p;
            self.occupant[p] = Some(q);
        }
    }

    /// Teleport the qubits whose core changes in `assignment`, before
    /// timeslice `t`. Each arriving qubit takes the free slot of its new core
    /// closest to its next partner.
    fn teleport(&mut self, assignment: &[usize], slices: &Timeslices, t: usize) {
        let k = self.k;
        let core_of = |p: usize| p / k;
        let arriving: Vec<usize> = (0..assignment.len())
            .filter(|&q| core_of(self.slot_of[q]) != assignment[q])
            .collect();
        // Free slots: empty, or left by a qubit that moves away.
        let mut free = vec![false; self.occupant.len()];
        for (p, f) in free.iter_mut().enumerate() {
            *f = self.occupant[p].is_none_or(|q| assignment[q] != core_of(p));
        }
        // Where each qubit will be: stayers now, arrivals once placed.
        let mut placed: Vec<Option<usize>> = (0..assignment.len())
            .map(|q| (core_of(self.slot_of[q]) == assignment[q]).then_some(self.slot_of[q]))
            .collect();
        let next_partner = |q: usize| (t..slices.len()).find_map(|s| slices.partner(s, q));
        let mut pending = arriving.clone();
        let mut moves = Vec::new();
        while !pending.is_empty() {
            // Arrivals whose partner is already placed in their core go first.
            let i = pending
                .iter()
                .position(|&q| {
                    next_partner(q)
                        .and_then(|p| placed[p])
                        .is_some_and(|p| core_of(p) == assignment[q])
                })
                .unwrap_or(0);
            let q = pending.remove(i);
            let core = assignment[q];
            let partner_slot = next_partner(q)
                .and_then(|p| placed[p])
                .filter(|&p| core_of(p) == core);
            let candidates = (core * k..(core + 1) * k).filter(|&p| free[p]);
            let local = |p: usize| p - core * k;
            let p = match partner_slot {
                _ if !self.options.partner_slots => candidates.min(),
                Some(s) => candidates.min_by_key(|&p| self.chip.distance(local(p), local(s))),
                // No partner there yet: the free slot with most free
                // neighbours, so that the partner can land next to it.
                None => candidates.max_by_key(|&p| {
                    let free_neighbours = self
                        .chip
                        .neighbours(local(p))
                        .iter()
                        .filter(|&&n| free[core * k + n])
                        .count();
                    (free_neighbours, usize::MAX - p)
                }),
            }
            .expect("the assignment fits the cores");
            free[p] = false;
            placed[q] = Some(p);
            moves.push((q, self.slot_of[q], p));
        }
        self.apply_moves(moves);
    }

    /// Emit teleportations `(qubit, from, to)`, all before the same timeslice,
    /// each when its destination is free if possible.
    fn apply_moves(&mut self, mut moves: Vec<(usize, usize, usize)>) {
        let k = self.k;
        while !moves.is_empty() {
            let i = moves
                .iter()
                .position(|&(_, _, to)| self.occupant[to].is_none())
                .unwrap_or(0);
            let (q, from, to) = moves.remove(i);
            self.out
                .push(Instruction::gate("teleport", vec![from, to], vec![]));
            self.teleportations += self.arch.hops(from / k, to / k);
            if self.occupant[from] == Some(q) {
                self.occupant[from] = None;
            }
            self.occupant[to] = Some(q);
            self.slot_of[q] = to;
        }
    }
}

/// The instructions of `work` on the qubits `members`, relabelled to
/// `0..members.len()`.
pub(super) fn local_circuit(work: &[&Instruction], members: &[usize]) -> Vec<Instruction> {
    let mut index = std::collections::HashMap::new();
    for (i, &q) in members.iter().enumerate() {
        index.insert(q, i);
    }
    work.iter()
        .filter(|inst| inst.qubits.iter().all(|q| index.contains_key(q)))
        .map(|inst| Instruction {
            qubits: inst.qubits.iter().map(|q| index[q]).collect(),
            ..(*inst).clone()
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gates::{matrix_1q, matrix_2q};
    use crate::linalg::{C, ZERO};
    use crate::multicore::{Method, Ports, Timing};
    use crate::testing::Rng as TestRng;

    /// Simulate `instructions` on `n` qubits from `state` (qubit 0 most
    /// significant). A run of teleportations moves states all at once.
    fn simulate(instructions: &[Instruction], n: usize, mut state: Vec<C>) -> Vec<C> {
        let bit = |q: usize| 1usize << (n - 1 - q);
        let mut i = 0;
        while i < instructions.len() {
            let inst = &instructions[i];
            if inst.name == "teleport" {
                let mut block = Vec::new();
                while i < instructions.len() && instructions[i].name == "teleport" {
                    block.push((instructions[i].qubits[0], instructions[i].qubits[1]));
                    i += 1;
                }
                let mut next = vec![ZERO; state.len()];
                for (index, &amplitude) in state.iter().enumerate() {
                    if amplitude == ZERO {
                        continue;
                    }
                    let mut moved = index;
                    for &(from, _) in &block {
                        moved &= !bit(from);
                    }
                    for &(from, to) in &block {
                        // Destinations not left by another qubit must be empty.
                        assert!(
                            block.iter().any(|&(f, _)| f == to) || index & bit(to) == 0,
                            "teleport into an occupied qubit {to}"
                        );
                        if index & bit(from) != 0 {
                            moved |= bit(to);
                        }
                    }
                    next[moved] += amplitude;
                }
                state = next;
                continue;
            }
            if inst.name == "barrier" {
                i += 1;
                continue;
            }
            let k = inst.qubits.len();
            let g: Vec<Vec<C>> = if k == 1 {
                matrix_1q(&inst.name, &inst.params)
                    .unwrap()
                    .0
                    .iter()
                    .map(|r| r.to_vec())
                    .collect()
            } else {
                matrix_2q(&inst.name, &inst.params)
                    .unwrap()
                    .0
                    .iter()
                    .map(|r| r.to_vec())
                    .collect()
            };
            let mut next = vec![ZERO; state.len()];
            for (index, &amplitude) in state.iter().enumerate() {
                if amplitude == ZERO {
                    continue;
                }
                let col = inst
                    .qubits
                    .iter()
                    .fold(0, |acc, &q| (acc << 1) | usize::from(index & bit(q) != 0));
                for (row, g_row) in g.iter().enumerate() {
                    let mut target = index;
                    for (pos, &q) in inst.qubits.iter().enumerate() {
                        if (row >> (k - 1 - pos)) & 1 == 1 {
                            target |= bit(q);
                        } else {
                            target &= !bit(q);
                        }
                    }
                    next[target] += g_row[col] * amplitude;
                }
            }
            state = next;
            i += 1;
        }
        state
    }

    /// The logical state of `qubits` placed by `layout` on `n` physical qubits.
    fn embed(logical: &[C], qubits: usize, layout: &[usize], n: usize) -> Vec<C> {
        let mut out = vec![ZERO; 1 << n];
        for (index, &amplitude) in logical.iter().enumerate() {
            let mut p = 0;
            for (q, &slot) in layout.iter().enumerate() {
                if index >> (qubits - 1 - q) & 1 == 1 {
                    p |= 1 << (n - 1 - slot);
                }
            }
            out[p] = amplitude;
        }
        out
    }

    fn random_circuit(rng: &mut TestRng, n: usize, gates: usize) -> Vec<Instruction> {
        let mut out: Vec<Instruction> = (0..n)
            .map(|q| {
                Instruction::gate(
                    "u",
                    vec![q],
                    vec![rng.range(0.0, 3.0), rng.range(0.0, 3.0), 0.4],
                )
            })
            .collect();
        for _ in 0..gates {
            if rng.below(3) == 0 {
                out.push(Instruction::gate(
                    "ry",
                    vec![rng.below(n)],
                    vec![rng.range(-3.0, 3.0)],
                ));
                continue;
            }
            let a = rng.below(n);
            let b = (a + 1 + rng.below(n - 1)) % n;
            out.push(match rng.below(3) {
                0 => Instruction::gate("cx", vec![a, b], vec![]),
                1 => Instruction::gate("cz", vec![a, b], vec![]),
                _ => Instruction::gate("rzz", vec![a, b], vec![0.7]),
            });
        }
        out
    }

    /// Default options with fewer layout trials, for speed.
    fn quick() -> PhysicalOptions {
        PhysicalOptions {
            sabre: SabreOptions {
                trials: 2,
                ..SabreOptions::default()
            },
            ..PhysicalOptions::default()
        }
    }

    /// Run the physical circuit and the original on the same input and
    /// compare, checking that gates and SWAPs stay on coupled qubits of one
    /// core and teleportations go between cores.
    fn check(
        circuit: &[Instruction],
        qubits: usize,
        arch: &Architecture,
        chip: &CouplingGraph,
        method: Method,
        options: &PhysicalOptions,
    ) {
        let (mapping, physical) = map_and_route(
            circuit,
            qubits,
            arch,
            chip,
            method,
            &Timing::default(),
            options,
        )
        .unwrap();
        let k = arch.qubits_per_core;
        for inst in &physical.instructions {
            if let [a, b] = inst.qubits[..] {
                if inst.name == "teleport" {
                    assert_ne!(a / k, b / k);
                    if let Some(ports) = &options.ports {
                        // From one port of a link to its other port.
                        assert_eq!(arch.hops(a / k, b / k), 1);
                        assert_eq!(a % k, ports.slot(&mapping.architecture, a / k, b / k));
                        assert_eq!(b % k, ports.slot(&mapping.architecture, b / k, a / k));
                    }
                } else {
                    assert_eq!(a / k, b / k, "{inst:?}");
                    assert!(chip.coupled(a % k, b % k), "{inst:?}");
                }
            }
        }
        let count = |name: &str| {
            physical
                .instructions
                .iter()
                .filter(|i| i.name == name)
                .count()
        };
        assert_eq!(count("swap"), physical.swaps);
        // Through ports, a move between unlinked cores is one teleportation
        // per link; otherwise it is one instruction.
        match options.ports {
            // Exchanges between full cores may add teleportations.
            Some(_) => {
                assert_eq!(count("teleport"), physical.teleportations);
                assert!(physical.teleportations >= mapping.communications());
            }
            None => {
                assert_eq!(count("teleport"), mapping.moved_qubits());
                assert_eq!(physical.teleportations, mapping.communications());
            }
        }
        assert_eq!(
            physical.instructions.len() - physical.swaps - count("teleport") - count("barrier"),
            circuit.iter().filter(|i| i.name != "barrier").count()
        );

        let n = arch.capacity();
        let mut zero = vec![ZERO; 1 << qubits];
        zero[0] = C::new(1.0, 0.0);
        let expected = simulate(circuit, qubits, zero);
        let start = embed(&[C::new(1.0, 0.0)], 0, &[], n);
        let got = simulate(&physical.instructions, n, start);
        let want = embed(&expected, qubits, &physical.final_layout, n);
        let overlap: C = got.iter().zip(&want).map(|(g, w)| w.conj() * g).sum();
        assert!(
            (overlap.norm() - 1.0).abs() < 1e-9,
            "overlap {}",
            overlap.norm()
        );
    }

    #[test]
    fn physical_circuits_are_equivalent_and_respect_the_chips() {
        let mut rng = TestRng::new(5);
        let cases = [
            (Architecture::new(3, 3).unwrap(), CouplingGraph::line(3)),
            (Architecture::line(3, 4).unwrap(), CouplingGraph::line(4)),
            (Architecture::new(2, 6).unwrap(), CouplingGraph::grid(2, 3)),
            (
                Architecture::ring(4, 3).unwrap(),
                CouplingGraph::all_to_all(3),
            ),
        ];
        for (arch, chip) in &cases {
            for trial in 0..6 {
                let qubits = arch.capacity() - trial % 3;
                let circuit = random_circuit(&mut rng, qubits, 40);
                for method in [Method::Hqa, Method::FgpRoee, Method::Naive] {
                    check(&circuit, qubits, arch, chip, method, &quick());
                }
                // Without the slot choice and the look-ahead, too.
                let plain = PhysicalOptions {
                    partner_slots: false,
                    look_ahead: 0,
                    ..quick()
                };
                check(&circuit, qubits, arch, chip, Method::Hqa, &plain);
            }
        }
    }

    #[test]
    fn circuits_through_ports_are_equivalent() {
        let mut rng = TestRng::new(9);
        let cases = [
            (Architecture::new(3, 4).unwrap(), CouplingGraph::line(4)),
            (Architecture::line(3, 5).unwrap(), CouplingGraph::line(5)),
            (Architecture::ring(4, 3).unwrap(), CouplingGraph::line(3)),
            (Architecture::grid(2, 2, 3).unwrap(), CouplingGraph::ring(3)),
        ];
        for (arch, chip) in &cases {
            let options = PhysicalOptions {
                ports: Some(Ports::central(arch, chip)),
                ..quick()
            };
            // The mapper keeps a slot free in every core.
            let most = arch.cores * (arch.qubits_per_core - 1);
            for trial in 0..5 {
                let qubits = most - trial % 2;
                let circuit = random_circuit(&mut rng, qubits, 40);
                for method in [
                    Method::Hqa,
                    Method::HqaLatency,
                    Method::FgpRoee,
                    Method::Naive,
                ] {
                    check(&circuit, qubits, arch, chip, method, &options);
                }
            }
        }
    }

    #[test]
    fn hqa_latency_runs_no_slower_than_hqa() {
        let mut rng = TestRng::new(8);
        let options = quick();
        for (arch, capacity) in [
            (Architecture::new(4, 4).unwrap(), 1),
            (Architecture::line(4, 4).unwrap(), 3),
        ] {
            let timing = Timing {
                link_capacity: capacity,
                ..Timing::default()
            };
            let chip = CouplingGraph::grid(2, 2);
            for _ in 0..4 {
                let circuit = random_circuit(&mut rng, 16, 120);
                let duration = |method| {
                    let (_, physical) =
                        map_and_route(&circuit, 16, &arch, &chip, method, &timing, &options)
                            .unwrap();
                    crate::multicore::schedule(&physical.instructions, &arch, &timing)
                        .unwrap()
                        .duration
                };
                assert!(duration(Method::HqaLatency) <= duration(Method::Hqa) + 1e-9);
            }
        }
    }

    #[test]
    fn single_qubit_circuits_and_measurements() {
        let arch = Architecture::new(2, 2).unwrap();
        let chip = CouplingGraph::line(2);
        let circuit = vec![
            Instruction::gate("h", vec![0], vec![]),
            Instruction::gate("barrier", vec![0, 1, 2], vec![]),
            Instruction::gate("x", vec![2], vec![]),
        ];
        check(&circuit, 3, &arch, &chip, Method::Hqa, &quick());
        let measured = vec![
            Instruction::gate("cx", vec![0, 1], vec![]),
            Instruction {
                name: "measure".into(),
                qubits: vec![1],
                params: vec![],
                clbits: vec![0],
            },
        ];
        let options = PhysicalOptions::default();
        let (_, physical) = map_and_route(
            &measured,
            2,
            &arch,
            &chip,
            Method::Hqa,
            &Timing::default(),
            &options,
        )
        .unwrap();
        assert_eq!(physical.instructions.last().unwrap().name, "measure");
    }

    #[test]
    fn errors() {
        let arch = Architecture::new(2, 3).unwrap();
        let options = PhysicalOptions::default();
        let wide = vec![Instruction::gate("ccx", vec![0, 1, 2], vec![])];
        assert!(
            map_and_route(
                &wide,
                3,
                &arch,
                &CouplingGraph::line(3),
                Method::Hqa,
                &Timing::default(),
                &options
            )
            .unwrap_err()
            .contains("decompose")
        );
        assert!(
            map_and_route(
                &[],
                3,
                &arch,
                &CouplingGraph::line(4),
                Method::Hqa,
                &Timing::default(),
                &options
            )
            .unwrap_err()
            .contains("the chip has 4 qubits")
        );
    }
}
