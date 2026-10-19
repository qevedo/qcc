//! SABRE: SWAP-based bidirectional heuristic search, from G. Li, Y. Ding and
//! Y. Xie, "Tackling the Qubit Mapping Problem for NISQ-Era Quantum Devices"
//! (ASPLOS 2019), with the release valve of Qiskit's implementation.
//!
//! Routing walks the circuit's dependency graph. Gates whose qubits are
//! coupled, and every single-qubit gate, run as soon as their predecessors
//! have. When every gate in the front layer is blocked, one SWAP is inserted:
//! the one, among SWAPs that touch a qubit of the front layer, that most
//! reduces
//!
//! ```text
//! H = max(decay(p₁), decay(p₂)) · ( Σ_front d / |front| + W · Σ_ext d / |ext| )
//! ```
//!
//! where `d` is the coupling distance between a gate's two qubits, `ext` is
//! the next gates after the front layer (the extended set), and `decay` grows
//! for qubits just swapped, so SWAPs spread out and can run in parallel. If
//! no gate runs for many SWAPs in a row (the release valve), the closest front
//! gate is routed along a shortest path.
//!
//! Not every SWAP costs three CNOTs. One right after a block of gates on the
//! same pair merges into it when the block is resynthesized: after a CX it
//! costs one more CNOT, after a generic three-CNOT block nothing. As in NASSC
//! (Liu, Li and Zhou, "Not All SWAPs Have the Same Cost", HPCA 2022), a SWAP
//! that brings the front layer closer has its score lowered by
//! `merge_weight` times the CNOTs it saves that way, and routings are
//! compared by the CNOTs their SWAPs cost.
//!
//! The initial layout is found by routing the circuit forwards and backwards
//! from a random layout, a few times: each pass starts from where the previous
//! one ended, so the layout settles where the start of the circuit wants it.
//! The best of several random starts is kept.

use super::coupling::CouplingGraph;
use crate::gates::{matrix_1q, matrix_2q};
use crate::ir::Instruction;
use crate::linalg::{M2, M4, kron};
use crate::rng::Rng;
use crate::two_qubit::cx_count;

/// Tuning of SABRE; the defaults follow the paper and Qiskit.
#[derive(Clone, Copy, Debug)]
pub struct SabreOptions {
    /// Random initial layouts tried by [`layout_and_route`].
    pub trials: usize,
    /// Forward–backward passes per trial.
    pub iterations: usize,
    /// Routings of each trial's layout, with different tie-breaking.
    pub swap_trials: usize,
    /// Gates in the extended set.
    pub extended_set: usize,
    /// Weight `W` of the extended set.
    pub extended_weight: f64,
    /// How much a SWAP raises the decay of its qubits.
    pub decay: f64,
    /// How much each CNOT a SWAP saves by merging with the block before it
    /// lowers its score (0: plain SABRE).
    pub merge_weight: f64,
    /// Whether [`layout_and_route`] first looks for a layout that needs no
    /// SWAP (see `perfect.rs`).
    pub perfect_layout: bool,
    pub seed: u64,
}

impl Default for SabreOptions {
    fn default() -> SabreOptions {
        SabreOptions {
            trials: 8,
            iterations: 3,
            swap_trials: 4,
            extended_set: 20,
            extended_weight: 0.5,
            decay: 0.001,
            merge_weight: 0.0,
            perfect_layout: false,
            seed: 1,
        }
    }
}

/// A routed circuit on physical qubits.
#[derive(Clone, Debug)]
pub struct Routed {
    /// The circuit on physical qubits, with SWAPs inserted (named "swap").
    pub instructions: Vec<Instruction>,
    /// `initial_layout[q]`: the physical qubit of logical qubit `q` at the start.
    pub initial_layout: Vec<usize>,
    /// `final_layout[q]`: where logical qubit `q` ends up.
    pub final_layout: Vec<usize>,
    pub swaps: usize,
    /// The CNOTs the SWAPs cost, merged with the blocks before them.
    pub swap_cost: usize,
}

/// Route `circuit` (on logical qubits `0..logical`) onto `graph`, starting
/// from `initial` (`initial[q]` is the physical qubit of logical qubit `q`).
pub fn route(
    circuit: &[Instruction],
    logical: usize,
    graph: &CouplingGraph,
    initial: &[usize],
    options: &SabreOptions,
) -> Result<Routed, String> {
    check(circuit, logical, graph)?;
    if initial.len() != logical {
        return Err(format!(
            "the layout places {} qubits, the circuit has {logical}",
            initial.len()
        ));
    }
    Router::new(circuit, logical, graph, initial, options).run()
}

/// [`route`], also looking at `ahead`: two-qubit gates (pairs of logical
/// qubits) that will run after `circuit`. They are not routed, but they join
/// the extended set when the circuit's own gates run out, so the last SWAPs
/// also prepare for what comes next.
pub fn route_ahead(
    circuit: &[Instruction],
    logical: usize,
    graph: &CouplingGraph,
    initial: &[usize],
    ahead: &[(usize, usize)],
    options: &SabreOptions,
) -> Result<Routed, String> {
    check(circuit, logical, graph)?;
    if initial.len() != logical || ahead.iter().any(|&(a, b)| a.max(b) >= logical) {
        return Err(format!(
            "the layout or the gates ahead do not match the {logical} qubits of the circuit"
        ));
    }
    let mut router = Router::new(circuit, logical, graph, initial, options);
    router.ahead = ahead;
    router.run()
}

/// Find an initial layout with SABRE's forward–backward search and route.
pub fn layout_and_route(
    circuit: &[Instruction],
    logical: usize,
    graph: &CouplingGraph,
    options: &SabreOptions,
) -> Result<Routed, String> {
    check(circuit, logical, graph)?;
    // Only two-qubit gates matter for the layout.
    let two_qubit: Vec<Instruction> = circuit
        .iter()
        .filter(|i| needs_coupling(i))
        .cloned()
        .collect();
    // A layout under which every gate's pair is coupled needs no SWAP.
    if options.perfect_layout {
        let pairs: Vec<(usize, usize)> = two_qubit
            .iter()
            .map(|i| (i.qubits[0], i.qubits[1]))
            .collect();
        if let Some(layout) = super::perfect::perfect_layout(&pairs, logical, graph) {
            return Router::new(circuit, logical, graph, &layout, options).run();
        }
    }
    let reversed: Vec<Instruction> = two_qubit.iter().rev().cloned().collect();
    let mut best: Option<Routed> = None;
    for trial in 0..options.trials.max(1) {
        let mut rng = Rng::new(options.seed.wrapping_add(trial as u64 * 0x9E37_79B9));
        let mut physical: Vec<usize> = (0..graph.qubits()).collect();
        for i in (1..physical.len()).rev() {
            physical.swap(i, rng.below(i + 1));
        }
        let mut layout = physical[..logical].to_vec();
        let pass_options = SabreOptions {
            seed: rng.next_u64(),
            ..*options
        };
        for _ in 0..options.iterations {
            layout = Router::new(&two_qubit, logical, graph, &layout, &pass_options)
                .run()?
                .final_layout;
            layout = Router::new(&reversed, logical, graph, &layout, &pass_options)
                .run()?
                .final_layout;
        }
        // Route the layout a few times with different tie-breaking, keep the best.
        for swap_trial in 0..options.swap_trials.max(1) {
            let route_options = SabreOptions {
                seed: pass_options.seed.wrapping_add(swap_trial as u64),
                ..*options
            };
            let routed = Router::new(circuit, logical, graph, &layout, &route_options).run()?;
            if best
                .as_ref()
                .is_none_or(|b| (routed.swap_cost, routed.swaps) < (b.swap_cost, b.swaps))
            {
                best = Some(routed);
            }
        }
    }
    Ok(best.expect("at least one trial"))
}

fn needs_coupling(inst: &Instruction) -> bool {
    inst.is_gate() && inst.qubits.len() == 2
}

fn check(circuit: &[Instruction], logical: usize, graph: &CouplingGraph) -> Result<(), String> {
    if logical > graph.qubits() {
        return Err(format!(
            "{logical} qubits do not fit on a chip of {}",
            graph.qubits()
        ));
    }
    for inst in circuit {
        if let Some(&q) = inst.qubits.iter().find(|&&q| q >= logical) {
            return Err(format!(
                "'{}' uses qubit {q}, beyond the {logical} of the circuit",
                inst.name
            ));
        }
        if inst.is_gate() && inst.qubits.len() > 2 {
            return Err(format!(
                "'{}' acts on {} qubits: decompose the circuit into one- and two-qubit gates first",
                inst.name,
                inst.qubits.len()
            ));
        }
    }
    Ok(())
}

struct Router<'a> {
    circuit: &'a [Instruction],
    graph: &'a CouplingGraph,
    options: &'a SabreOptions,
    /// `physical[q]`: where logical qubit `q` is now.
    physical: Vec<usize>,
    /// `logical[p]`: the logical qubit on physical qubit `p`, if any.
    logical: Vec<Option<usize>>,
    decay: Vec<f64>,
    rng: Rng,
    out: Vec<Instruction>,
    swaps: usize,
    initial: Vec<usize>,
    /// Gates after the circuit, for the extended set only.
    ahead: &'a [(usize, usize)],
    /// The open block of gates on each physical qubit's pair, if any.
    open: Vec<Option<usize>>,
    blocks: Vec<OpenBlock>,
    swap_cost: usize,
}

/// Gates on one pair of physical qubits since either last met another qubit.
struct OpenBlock {
    pair: [usize; 2],
    /// The block's unitary, `pair[0]` most significant; None once a gate has
    /// no known matrix.
    matrix: Option<M4>,
    /// CNOTs saved by a SWAP right after it, once computed.
    saving: Option<usize>,
}

impl<'a> Router<'a> {
    fn new(
        circuit: &'a [Instruction],
        logical: usize,
        graph: &'a CouplingGraph,
        initial: &[usize],
        options: &'a SabreOptions,
    ) -> Router<'a> {
        let mut at = vec![None; graph.qubits()];
        for (q, &p) in initial.iter().enumerate() {
            at[p] = Some(q);
        }
        debug_assert!(initial.len() == logical);
        Router {
            circuit,
            graph,
            options,
            physical: initial.to_vec(),
            logical: at,
            decay: vec![1.0; graph.qubits()],
            rng: Rng::new(options.seed),
            out: Vec::with_capacity(circuit.len()),
            swaps: 0,
            initial: initial.to_vec(),
            ahead: &[],
            open: vec![None; graph.qubits()],
            blocks: Vec::new(),
            swap_cost: 0,
        }
    }

    fn run(mut self) -> Result<Routed, String> {
        let n = self.circuit.len();
        // The dependency graph: each instruction follows the previous one on
        // each of its qubits.
        let mut remaining = vec![0usize; n];
        let mut successors: Vec<Vec<usize>> = vec![Vec::new(); n];
        let mut last: Vec<Option<usize>> = vec![None; self.physical.len()];
        for (i, inst) in self.circuit.iter().enumerate() {
            for &q in &inst.qubits {
                if let Some(j) = last[q] {
                    if !successors[j].contains(&i) {
                        successors[j].push(i);
                        remaining[i] += 1;
                    }
                }
                last[q] = Some(i);
            }
        }
        let mut front: Vec<usize> = (0..n).filter(|&i| remaining[i] == 0).collect();
        let mut done = 0;
        let mut stalled = 0;
        let release_after = 10 * self.graph.qubits().max(10);
        loop {
            // Run everything that can run.
            let mut progressed = false;
            loop {
                let mut next = Vec::with_capacity(front.len());
                let mut ran = false;
                for &node in &front {
                    if self.can_run(node) {
                        self.emit(node);
                        done += 1;
                        ran = true;
                        for &s in &successors[node] {
                            remaining[s] -= 1;
                            if remaining[s] == 0 {
                                next.push(s);
                            }
                        }
                    } else {
                        next.push(node);
                    }
                }
                front = next;
                progressed |= ran;
                if !ran {
                    break;
                }
            }
            if done == n {
                break;
            }
            if progressed {
                self.decay.iter_mut().for_each(|d| *d = 1.0);
                stalled = 0;
            }
            if stalled >= release_after {
                self.release_valve(&front);
                stalled = 0;
                continue;
            }
            let extended = self.extended_set(&front, &successors, &remaining);
            let (a, b) = self.best_swap(&front, &extended);
            self.swap(a, b);
            stalled += 1;
            if self.swaps.is_multiple_of(5) {
                self.decay.iter_mut().for_each(|d| *d = 1.0);
            }
        }
        Ok(Routed {
            instructions: self.out,
            initial_layout: self.initial,
            final_layout: self.physical,
            swaps: self.swaps,
            swap_cost: self.swap_cost,
        })
    }

    fn can_run(&self, node: usize) -> bool {
        let inst = &self.circuit[node];
        !needs_coupling(inst)
            || self
                .graph
                .coupled(self.physical[inst.qubits[0]], self.physical[inst.qubits[1]])
    }

    fn emit(&mut self, node: usize) {
        let inst = &self.circuit[node];
        let out = Instruction {
            name: inst.name.clone(),
            qubits: inst.qubits.iter().map(|&q| self.physical[q]).collect(),
            params: inst.params.clone(),
            clbits: inst.clbits.clone(),
        };
        if self.options.merge_weight > 0.0 {
            self.track(&out);
        }
        self.out.push(out);
    }

    /// Add `inst` (on physical qubits) to the open blocks.
    fn track(&mut self, inst: &Instruction) {
        match *inst.qubits.as_slice() {
            [p] if inst.is_gate() => {
                if let Some(b) = self.open[p] {
                    let block = &mut self.blocks[b];
                    block.matrix = match (block.matrix, matrix_1q(&inst.name, &inst.params)) {
                        (Some(m), Some(g)) if block.pair[0] == p => {
                            Some(kron(&g, &M2::IDENTITY) * m)
                        }
                        (Some(m), Some(g)) => Some(kron(&M2::IDENTITY, &g) * m),
                        _ => None,
                    };
                    block.saving = None;
                }
            }
            [p, q] if inst.is_gate() => {
                let g = matrix_2q(&inst.name, &inst.params);
                match (self.open[p], self.open[q]) {
                    (Some(b), Some(c)) if b == c => {
                        let block = &mut self.blocks[b];
                        let oriented = g.map(|g| orient(g, block.pair[0] == p));
                        block.matrix = match (block.matrix, oriented) {
                            (Some(m), Some(g)) => Some(g * m),
                            _ => None,
                        };
                        block.saving = None;
                    }
                    _ => {
                        self.blocks.push(OpenBlock {
                            pair: [p, q],
                            matrix: g,
                            saving: None,
                        });
                        self.open[p] = Some(self.blocks.len() - 1);
                        self.open[q] = Some(self.blocks.len() - 1);
                    }
                }
            }
            _ => {
                for &p in &inst.qubits {
                    self.open[p] = None;
                }
            }
        }
    }

    /// CNOTs a SWAP on `(a, b)` would save by merging with the block before it.
    fn saving(&mut self, a: usize, b: usize) -> usize {
        let (Some(x), Some(y)) = (self.open[a], self.open[b]) else {
            return 0;
        };
        if x != y {
            return 0;
        }
        let block = &mut self.blocks[x];
        if block.saving.is_none() {
            block.saving = Some(match block.matrix {
                Some(m) => {
                    let swap = matrix_2q("swap", &[]).expect("a standard gate");
                    (cx_count(&m) + 3).saturating_sub(cx_count(&(swap * m)))
                }
                None => 0,
            });
        }
        block.saving.unwrap()
    }

    fn swap(&mut self, a: usize, b: usize) {
        let swap = Instruction::gate("swap", vec![a, b], vec![]);
        if self.options.merge_weight > 0.0 {
            self.swap_cost += 3 - self.saving(a, b).min(3);
            self.track(&swap);
        } else {
            self.swap_cost += 3;
        }
        self.out.push(swap);
        let (qa, qb) = (self.logical[a], self.logical[b]);
        self.logical[a] = qb;
        self.logical[b] = qa;
        if let Some(q) = qa {
            self.physical[q] = b;
        }
        if let Some(q) = qb {
            self.physical[q] = a;
        }
        self.decay[a] += self.options.decay;
        self.decay[b] += self.options.decay;
        self.swaps += 1;
    }

    /// The qubit pairs of the next two-qubit gates after the front layer.
    fn extended_set(
        &self,
        front: &[usize],
        successors: &[Vec<usize>],
        remaining: &[usize],
    ) -> Vec<(usize, usize)> {
        let mut out = Vec::new();
        let mut seen = vec![false; self.circuit.len()];
        let mut queue: std::collections::VecDeque<usize> = front.iter().copied().collect();
        // Track the predecessors met so far, so only gates whose blocking
        // predecessors are all in the front layer or the set count early.
        let mut pending = remaining.to_vec();
        while let Some(node) = queue.pop_front() {
            for &s in &successors[node] {
                pending[s] = pending[s].saturating_sub(1);
                if pending[s] == 0 && !seen[s] {
                    seen[s] = true;
                    let inst = &self.circuit[s];
                    if needs_coupling(inst) {
                        out.push((inst.qubits[0], inst.qubits[1]));
                        if out.len() >= self.options.extended_set {
                            return out;
                        }
                    }
                    queue.push_back(s);
                }
            }
        }
        let room = self.options.extended_set - out.len();
        out.extend(self.ahead.iter().take(room));
        out
    }

    fn best_swap(&mut self, front: &[usize], extended: &[(usize, usize)]) -> (usize, usize) {
        let gates: Vec<(usize, usize)> = front
            .iter()
            .map(|&g| &self.circuit[g])
            .filter(|i| needs_coupling(i))
            .map(|i| (i.qubits[0], i.qubits[1]))
            .collect();
        let mut candidates: Vec<(usize, usize)> = Vec::new();
        for &(a, b) in &gates {
            for q in [a, b] {
                let p = self.physical[q];
                for &nb in self.graph.neighbours(p) {
                    let edge = (p.min(nb), p.max(nb));
                    if !candidates.contains(&edge) {
                        candidates.push(edge);
                    }
                }
            }
        }
        let mut best: Vec<(usize, usize)> = Vec::new();
        let mut best_score = f64::INFINITY;
        let front_now: f64 = gates
            .iter()
            .map(|&(a, b)| self.graph.distance(self.physical[a], self.physical[b]) as f64)
            .sum();
        for &(p1, p2) in &candidates {
            let at = |q: usize| {
                let p = self.physical[q];
                if p == p1 {
                    p2
                } else if p == p2 {
                    p1
                } else {
                    p
                }
            };
            let cost = |pairs: &[(usize, usize)]| -> f64 {
                if pairs.is_empty() {
                    return 0.0;
                }
                pairs
                    .iter()
                    .map(|&(a, b)| self.graph.distance(at(a), at(b)) as f64)
                    .sum::<f64>()
                    / pairs.len() as f64
            };
            let mut score = self.decay[p1].max(self.decay[p2])
                * (cost(&gates) + self.options.extended_weight * cost(extended));
            // Merging is a bonus only for SWAPs that bring the front closer:
            // a free SWAP that does not would be chosen again and again.
            let closer = cost(&gates) * (gates.len() as f64) < front_now - 1e-9;
            if self.options.merge_weight > 0.0 && closer {
                score -= self.options.merge_weight * self.saving(p1, p2) as f64;
            }
            if score < best_score - 1e-12 {
                best_score = score;
                best.clear();
                best.push((p1, p2));
            } else if (score - best_score).abs() <= 1e-12 {
                best.push((p1, p2));
            }
        }
        best[self.rng.below(best.len())]
    }

    /// Route the closest front gate along a shortest path.
    fn release_valve(&mut self, front: &[usize]) {
        let (a, b) = front
            .iter()
            .map(|&g| &self.circuit[g])
            .filter(|i| needs_coupling(i))
            .map(|i| (self.physical[i.qubits[0]], self.physical[i.qubits[1]]))
            .min_by_key(|&(a, b)| self.graph.distance(a, b))
            .expect("a blocked gate");
        let path = self.graph.path(a, b);
        for w in path.windows(2).take(path.len().saturating_sub(2)) {
            self.swap(w[0], w[1]);
        }
    }
}

/// Gate `g` (on its operands in order) in a block whose first qubit is its
/// first operand (`forward`) or its second.
fn orient(g: M4, forward: bool) -> M4 {
    if forward {
        g
    } else {
        let swap = matrix_2q("swap", &[]).expect("a standard gate");
        swap * g * swap
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testing::{Rng as TestRng, circuit_unitary, dense_equal_up_to_phase};

    fn random_circuit(rng: &mut TestRng, n: usize, gates: usize) -> Vec<Instruction> {
        (0..gates)
            .map(|_| {
                if rng.below(3) == 0 {
                    return Instruction::gate("rz", vec![rng.below(n)], vec![rng.range(-3.0, 3.0)]);
                }
                let a = rng.below(n);
                let b = (a + 1 + rng.below(n - 1)) % n;
                match rng.below(3) {
                    0 => Instruction::gate("cx", vec![a, b], vec![]),
                    1 => Instruction::gate("cz", vec![a, b], vec![]),
                    _ => Instruction::gate("rzz", vec![a, b], vec![0.7]),
                }
            })
            .collect()
    }

    /// The permutation that puts logical qubit `q` on physical qubit `layout[q]`,
    /// as a function on basis indices (qubit 0 most significant).
    fn permute(index: usize, layout: &[usize], n: usize) -> usize {
        let mut out = 0;
        for (q, &p) in layout.iter().enumerate() {
            if index >> (n - 1 - q) & 1 == 1 {
                out |= 1 << (n - 1 - p);
            }
        }
        out
    }

    /// U_routed · P_initial = P_final · U_logical.
    fn check_equivalent(logical: &[Instruction], routed: &Routed, n: usize) {
        let u = circuit_unitary(logical, n);
        let r = circuit_unitary(&routed.instructions, n);
        let dim = 1 << n;
        // Columns: r applied to P_initial |j⟩ must be P_final applied to u |j⟩.
        let mut lhs = vec![vec![crate::linalg::ZERO; dim]; dim];
        let mut rhs = vec![vec![crate::linalg::ZERO; dim]; dim];
        for j in 0..dim {
            let pj = permute(j, &routed.initial_layout, n);
            for i in 0..dim {
                lhs[i][j] = r[i][pj];
                rhs[permute(i, &routed.final_layout, n)][j] = u[i][j];
            }
        }
        assert!(dense_equal_up_to_phase(&lhs, &rhs, 1e-9));
    }

    #[test]
    fn routed_circuits_respect_the_coupling_and_are_equivalent() {
        let mut rng = TestRng::new(2);
        for graph in [
            CouplingGraph::line(5),
            CouplingGraph::ring(6),
            CouplingGraph::grid(2, 3),
            CouplingGraph::line(4),
        ] {
            let n = graph.qubits();
            for trial in 0..10 {
                let circuit = random_circuit(&mut rng, n, 30);
                let options = SabreOptions {
                    seed: trial,
                    ..SabreOptions::default()
                };
                for routed in [
                    route(&circuit, n, &graph, &(0..n).collect::<Vec<_>>(), &options).unwrap(),
                    layout_and_route(&circuit, n, &graph, &options).unwrap(),
                ] {
                    for inst in &routed.instructions {
                        if inst.qubits.len() == 2 {
                            assert!(graph.coupled(inst.qubits[0], inst.qubits[1]), "{inst:?}");
                        }
                    }
                    let non_swaps = routed
                        .instructions
                        .iter()
                        .filter(|i| i.name != "swap")
                        .count();
                    assert_eq!(
                        (non_swaps, routed.instructions.len() - non_swaps),
                        (circuit.len(), routed.swaps)
                    );
                    check_equivalent(&circuit, &routed, n);
                }
            }
        }
    }

    #[test]
    fn gates_ahead_steer_but_are_not_routed() {
        let mut rng = TestRng::new(6);
        let graph = CouplingGraph::grid(2, 3);
        let identity: Vec<usize> = (0..6).collect();
        for seed in 0..5 {
            let circuit = random_circuit(&mut rng, 6, 20);
            let ahead: Vec<(usize, usize)> = (0..6).map(|q| (q, (q + 3) % 6)).collect();
            let options = SabreOptions {
                seed,
                ..SabreOptions::default()
            };
            let routed = route_ahead(&circuit, 6, &graph, &identity, &ahead, &options).unwrap();
            assert_eq!(routed.instructions.len(), circuit.len() + routed.swaps);
            check_equivalent(&circuit, &routed, 6);
        }
        assert!(
            route_ahead(
                &[],
                6,
                &graph,
                &identity,
                &[(0, 6)],
                &SabreOptions::default()
            )
            .is_err()
        );
    }

    #[test]
    fn fewer_logical_than_physical_qubits() {
        let mut rng = TestRng::new(3);
        let graph = CouplingGraph::grid(3, 4);
        let circuit = random_circuit(&mut rng, 7, 200);
        let routed = layout_and_route(&circuit, 7, &graph, &SabreOptions::default()).unwrap();
        assert!(
            routed
                .instructions
                .iter()
                .all(|i| i.qubits.len() < 2 || graph.coupled(i.qubits[0], i.qubits[1]))
        );
        let mut used = routed.final_layout.clone();
        used.sort();
        used.dedup();
        assert_eq!(used.len(), 7);
    }

    #[test]
    fn the_layout_search_beats_a_fixed_layout_on_average() {
        let mut rng = TestRng::new(4);
        let graph = CouplingGraph::grid(4, 4);
        let (mut fixed, mut searched) = (0, 0);
        for seed in 0..6 {
            let circuit = random_circuit(&mut rng, 16, 150);
            let options = SabreOptions {
                seed,
                ..SabreOptions::default()
            };
            fixed += route(&circuit, 16, &graph, &(0..16).collect::<Vec<_>>(), &options)
                .unwrap()
                .swaps;
            searched += layout_and_route(&circuit, 16, &graph, &options)
                .unwrap()
                .swaps;
        }
        assert!(searched < fixed, "searched {searched}, fixed {fixed}");
    }

    #[test]
    fn errors() {
        let graph = CouplingGraph::line(3);
        let wide = vec![Instruction::gate("ccx", vec![0, 1, 2], vec![])];
        assert!(
            layout_and_route(&wide, 3, &graph, &SabreOptions::default())
                .unwrap_err()
                .contains("decompose")
        );
        assert!(
            layout_and_route(&[], 4, &graph, &SabreOptions::default())
                .unwrap_err()
                .contains("do not fit")
        );
    }
}
