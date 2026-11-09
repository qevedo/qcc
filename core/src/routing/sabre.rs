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
//!
//! One search routes the same circuit hundreds of times (32 starts, each
//! with six layout passes and 16 routings on the device), so routing is
//! written to cost little per SWAP while making exactly the choices it
//! always made: the dependency graph and every gate's matrix are worked out
//! once per circuit ([`Dag`]); the extended set reuses its scratch space
//! instead of allocating two arrays as long as the circuit for every SWAP;
//! a routing is recorded as the order of its nodes and SWAPs and written out
//! as instructions only if it is kept ([`Trace`]); and the CNOTs a SWAP
//! saves after a block are remembered by the block's unitary across the
//! routings of one start ([`Savings`]). On four 16-qubit benchmark circuits
//! (Quantum Volume, QFT, Clifford+T, Hamiltonian) one search takes 2.1–2.9×
//! fewer instructions for the same routings, and with the same options every
//! benchmark compilation gives the same output as before.

use super::coupling::CouplingGraph;
use crate::gates::{matrix_1q, matrix_2q};
use crate::ir::Instruction;
use crate::linalg::{M2, M4};
use crate::rng::Rng;
use crate::two_qubit::cx_count;
use std::cell::RefCell;
use std::collections::HashMap;

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
    let dag = Dag::new(circuit, logical, options);
    Ok(
        Router::new(&dag, &Savings::default(), graph, initial, options)
            .run()
            .into_routed(circuit),
    )
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
    let dag = Dag::new(circuit, logical, options);
    let savings = Savings::default();
    let mut router = Router::new(&dag, &savings, graph, initial, options);
    router.ahead = ahead;
    Ok(router.run().into_routed(circuit))
}

/// Find an initial layout with SABRE's forward–backward search and route.
pub fn layout_and_route(
    circuit: &[Instruction],
    logical: usize,
    graph: &CouplingGraph,
    options: &SabreOptions,
) -> Result<Routed, String> {
    let mut best = layout_and_route_and_first(circuit, logical, graph, options, 0)?;
    Ok(best.swap_remove(0))
}

/// [`layout_and_route`]'s routing, and also the one it would have chosen
/// with only the first `first` swap trials of each layout, if that differs
/// (and `first` is less than `options.swap_trials`): the best of all, then
/// that one. Device compilation optimizes both, as the CNOTs the SWAPs cost
/// are an estimate of what the optimized circuit will need, and keeps the
/// cheaper, so searching more widely can't make a circuit worse than the
/// narrower search would have.
pub fn layout_and_route_and_first(
    circuit: &[Instruction],
    logical: usize,
    graph: &CouplingGraph,
    options: &SabreOptions,
    first: usize,
) -> Result<Vec<Routed>, String> {
    check(circuit, logical, graph)?;
    // Only two-qubit gates matter for the layout.
    let two_qubit: Vec<Instruction> = circuit
        .iter()
        .filter(|i| needs_coupling(i))
        .cloned()
        .collect();
    // Every routing below of the same circuit shares its dependency graph.
    let whole = Dag::new(circuit, logical, options);
    // A layout under which every gate's pair is coupled needs no SWAP.
    if options.perfect_layout {
        let pairs: Vec<(usize, usize)> = two_qubit
            .iter()
            .map(|i| (i.qubits[0], i.qubits[1]))
            .collect();
        if let Some(layout) = super::perfect::perfect_layout(&pairs, logical, graph) {
            return Ok(vec![
                Router::new(&whole, &Savings::default(), graph, &layout, options)
                    .run()
                    .into_routed(circuit),
            ]);
        }
    }
    let reversed: Vec<Instruction> = two_qubit.iter().rev().cloned().collect();
    let forwards = Dag::new(&two_qubit, logical, options);
    let backwards = Dag::new(&reversed, logical, options);
    // Each trial is seeded by its number alone, so the trials run side by
    // side and come back in order; the best is then chosen as one loop over
    // them would (the first of the cheapest), whatever the threads did.
    let trial = |trial: usize| -> (Trace, Option<Trace>) {
        let mut rng = Rng::new(options.seed.wrapping_add(trial as u64 * 0x9E37_79B9));
        let mut physical: Vec<usize> = (0..graph.qubits()).collect();
        for i in (1..physical.len()).rev() {
            physical.swap(i, rng.below(i + 1));
        }
        let mut layout = physical[..logical].to_vec();
        // One trial runs on one thread: its routings share what they learn
        // about blocks.
        let savings = Savings::default();
        let pass_options = SabreOptions {
            seed: rng.next_u64(),
            ..*options
        };
        for _ in 0..options.iterations {
            layout = Router::new(&forwards, &savings, graph, &layout, &pass_options)
                .run()
                .final_layout;
            layout = Router::new(&backwards, &savings, graph, &layout, &pass_options)
                .run()
                .final_layout;
        }
        // Route the layout a few times with different tie-breaking, keep the
        // best, and the best of the first few as it stood after them.
        let mut best: Option<Trace> = None;
        let mut best_first: Option<Trace> = None;
        for swap_trial in 0..options.swap_trials.max(1) {
            if swap_trial == first {
                best_first = best.clone();
            }
            let route_options = SabreOptions {
                seed: pass_options.seed.wrapping_add(swap_trial as u64),
                ..*options
            };
            let routed = Router::new(&whole, &savings, graph, &layout, &route_options).run();
            if best.as_ref().is_none_or(|b| routed.better_than(b)) {
                best = Some(routed);
            }
        }
        (best.expect("at least one swap trial"), best_first)
    };
    let mut best: Option<Trace> = None;
    let mut best_first: Option<Trace> = None;
    for (routed, first_routed) in crate::parallel::map_ordered(options.trials.max(1), trial) {
        if let Some(f) = first_routed {
            if best_first.as_ref().is_none_or(|b| f.better_than(b)) {
                best_first = Some(f);
            }
        }
        if best.as_ref().is_none_or(|b| routed.better_than(b)) {
            best = Some(routed);
        }
    }
    // Only the routings kept are written out as instructions.
    let best = best.expect("at least one trial");
    let mut out = vec![];
    if let Some(f) = best_first.filter(|f| f.steps != best.steps || f.initial != best.initial) {
        out.push(f.into_routed(circuit));
    }
    out.insert(0, best.into_routed(circuit));
    Ok(out)
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

/// What a routing needs to know about a circuit, worked out once and shared
/// by all its routings: one search routes the same circuit hundreds of
/// times (layout passes and swap trials), and building the dependency graph
/// and every gate's matrix again for each was a large part of the time.
struct Dag<'a> {
    circuit: &'a [Instruction],
    /// Physical qubits are numbered beyond the logical ones; the dependency
    /// graph only needs the logical ones.
    logical: usize,
    /// The successors of node `i` are `successors[starts[i]..starts[i + 1]]`:
    /// each instruction follows the previous one on each of its qubits.
    starts: Vec<usize>,
    successors: Vec<usize>,
    /// How many predecessors each node waits for.
    predecessors: Vec<usize>,
    /// Each node's kind.
    kinds: Vec<Kind>,
    /// The matrices of the gates, for the open blocks (only when SWAPs that
    /// merge are scored): `Kind::OneQubit(Some(i))` is `ones[i]`, and
    /// `Kind::TwoQubit(Some(i))` is `twos[i]`, on its operands in order and
    /// then the other way round.
    ones: Vec<M2>,
    twos: Vec<[M4; 2]>,
}

#[derive(Clone, Copy)]
enum Kind {
    /// A single-qubit gate, with its matrix if it has a known one.
    OneQubit(Option<u32>),
    /// A two-qubit gate, with its matrices if known.
    TwoQubit(Option<u32>),
    /// Anything else: measurements, resets, barriers. It closes the blocks
    /// of the qubits it touches.
    Other,
}

impl<'a> Dag<'a> {
    fn new(circuit: &'a [Instruction], logical: usize, options: &SabreOptions) -> Dag<'a> {
        let n = circuit.len();
        let mut lists: Vec<Vec<usize>> = vec![Vec::new(); n];
        let mut predecessors = vec![0usize; n];
        let mut last: Vec<Option<usize>> = vec![None; logical];
        for (i, inst) in circuit.iter().enumerate() {
            for &q in &inst.qubits {
                if let Some(j) = last[q] {
                    if !lists[j].contains(&i) {
                        lists[j].push(i);
                        predecessors[i] += 1;
                    }
                }
                last[q] = Some(i);
            }
        }
        let mut starts = Vec::with_capacity(n + 1);
        let mut successors = Vec::with_capacity(lists.iter().map(Vec::len).sum());
        for list in &lists {
            starts.push(successors.len());
            successors.extend_from_slice(list);
        }
        starts.push(successors.len());
        let tracked = options.merge_weight > 0.0;
        let (mut ones, mut twos) = (Vec::new(), Vec::new());
        let kinds = circuit
            .iter()
            .map(|inst| match inst.qubits.len() {
                1 if inst.is_gate() => Kind::OneQubit(
                    tracked
                        .then(|| matrix_1q(&inst.name, &inst.params))
                        .flatten()
                        .map(|m| {
                            ones.push(m);
                            ones.len() as u32 - 1
                        }),
                ),
                2 if inst.is_gate() => Kind::TwoQubit(
                    tracked
                        .then(|| matrix_2q(&inst.name, &inst.params))
                        .flatten()
                        .map(|m| {
                            twos.push([m, reversed(&m)]);
                            twos.len() as u32 - 1
                        }),
                ),
                _ => Kind::Other,
            })
            .collect();
        Dag {
            circuit,
            logical,
            starts,
            successors,
            predecessors,
            kinds,
            ones,
            twos,
        }
    }

    fn successors(&self, node: usize) -> &[usize] {
        &self.successors[self.starts[node]..self.starts[node + 1]]
    }
}

/// One routing as it ran: the circuit's nodes and the SWAPs in output
/// order. It becomes instructions only if it is kept ([`Trace::into_routed`]):
/// the search writes hundreds of routings and keeps one.
#[derive(Clone)]
struct Trace {
    steps: Vec<Step>,
    initial: Vec<usize>,
    final_layout: Vec<usize>,
    swaps: usize,
    swap_cost: usize,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Step {
    /// Node `i` of the circuit, on wherever its qubits are at that point.
    Node(usize),
    /// A SWAP of two physical qubits.
    Swap(usize, usize),
}

impl Trace {
    /// Whether this routing costs less than `other`: fewer CNOTs for its
    /// SWAPs, then fewer SWAPs.
    fn better_than(&self, other: &Trace) -> bool {
        (self.swap_cost, self.swaps) < (other.swap_cost, other.swaps)
    }

    /// The routed circuit, replaying the SWAPs to place each node.
    fn into_routed(self, circuit: &[Instruction]) -> Routed {
        let mut physical = self.initial.clone();
        let size = physical.iter().map(|&p| p + 1).max().unwrap_or(0);
        let size = self
            .steps
            .iter()
            .map(|s| match *s {
                Step::Swap(a, b) => a.max(b) + 1,
                Step::Node(_) => 0,
            })
            .max()
            .unwrap_or(0)
            .max(size);
        let mut at: Vec<Option<usize>> = vec![None; size];
        for (q, &p) in physical.iter().enumerate() {
            at[p] = Some(q);
        }
        let mut instructions = Vec::with_capacity(self.steps.len());
        for step in &self.steps {
            match *step {
                Step::Node(node) => {
                    let inst = &circuit[node];
                    instructions.push(Instruction {
                        name: inst.name.clone(),
                        qubits: inst.qubits.iter().map(|&q| physical[q]).collect(),
                        params: inst.params.clone(),
                        clbits: inst.clbits.clone(),
                    });
                }
                Step::Swap(a, b) => {
                    instructions.push(Instruction::gate("swap", vec![a, b], vec![]));
                    let (qa, qb) = (at[a], at[b]);
                    at[a] = qb;
                    at[b] = qa;
                    if let Some(q) = qa {
                        physical[q] = b;
                    }
                    if let Some(q) = qb {
                        physical[q] = a;
                    }
                }
            }
        }
        debug_assert_eq!(physical, self.final_layout);
        Routed {
            instructions,
            initial_layout: self.initial,
            final_layout: self.final_layout,
            swaps: self.swaps,
            swap_cost: self.swap_cost,
        }
    }
}

/// The CNOTs a SWAP saves after a block, by the block's unitary, bit for
/// bit: the routings of one layout keep meeting the same blocks, and two
/// CNOT counts (two KAK decompositions) per block were a third of the
/// routing's time. A unitary equal in value but not in bits is only a miss.
type Savings = RefCell<HashMap<[u64; 32], usize>>;

fn bits(m: &M4) -> [u64; 32] {
    let mut out = [0; 32];
    for (k, c) in m.0.iter().flatten().enumerate() {
        out[2 * k] = c.re.to_bits();
        out[2 * k + 1] = c.im.to_bits();
    }
    out
}

struct Router<'a> {
    dag: &'a Dag<'a>,
    savings: &'a Savings,
    graph: &'a CouplingGraph,
    options: &'a SabreOptions,
    /// `physical[q]`: where logical qubit `q` is now.
    physical: Vec<usize>,
    /// `logical[p]`: the logical qubit on physical qubit `p`, if any.
    logical: Vec<Option<usize>>,
    decay: Vec<f64>,
    rng: Rng,
    steps: Vec<Step>,
    swaps: usize,
    initial: Vec<usize>,
    /// Gates after the circuit, for the extended set only.
    ahead: &'a [(usize, usize)],
    /// The open block of gates on each physical qubit's pair, if any.
    open: Vec<Option<usize>>,
    blocks: Vec<OpenBlock>,
    swap_cost: usize,
    /// Scratch space for [`Router::extended_set`], kept between calls so
    /// that each call touches only the nodes it visits: a node's entries in
    /// `pending` and `seen` count only if its stamp is the current call's.
    stamp: u32,
    pending_stamp: Vec<u32>,
    pending: Vec<usize>,
    seen: Vec<u32>,
    queue: std::collections::VecDeque<usize>,
    extended: Vec<(usize, usize)>,
    /// Scratch space for [`Router::best_swap`]: the front's and the
    /// extended set's pairs on physical qubits, and the SWAPs it weighs.
    gates: Vec<(usize, usize)>,
    extended_on_chip: Vec<(usize, usize)>,
    candidates: Vec<(usize, usize)>,
    best: Vec<(usize, usize)>,
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
        dag: &'a Dag<'a>,
        savings: &'a Savings,
        graph: &'a CouplingGraph,
        initial: &[usize],
        options: &'a SabreOptions,
    ) -> Router<'a> {
        let mut at = vec![None; graph.qubits()];
        for (q, &p) in initial.iter().enumerate() {
            at[p] = Some(q);
        }
        debug_assert!(initial.len() == dag.logical);
        let n = dag.circuit.len();
        Router {
            dag,
            savings,
            graph,
            options,
            physical: initial.to_vec(),
            logical: at,
            decay: vec![1.0; graph.qubits()],
            rng: Rng::new(options.seed),
            steps: Vec::with_capacity(n),
            swaps: 0,
            initial: initial.to_vec(),
            ahead: &[],
            open: vec![None; graph.qubits()],
            // A block opens at most once per two-qubit gate or SWAP: room
            // for the gates and as many SWAPs, so it rarely has to grow.
            blocks: Vec::with_capacity(if options.merge_weight > 0.0 {
                2 * dag.twos.len()
            } else {
                0
            }),
            swap_cost: 0,
            stamp: 0,
            pending_stamp: vec![0; n],
            pending: vec![0; n],
            seen: vec![0; n],
            queue: std::collections::VecDeque::new(),
            extended: Vec::new(),
            gates: Vec::new(),
            extended_on_chip: Vec::new(),
            candidates: Vec::new(),
            best: Vec::new(),
        }
    }

    fn run(mut self) -> Trace {
        let dag = self.dag;
        let n = dag.circuit.len();
        let mut remaining = dag.predecessors.clone();
        let mut front: Vec<usize> = (0..n).filter(|&i| remaining[i] == 0).collect();
        let mut next = Vec::with_capacity(front.len());
        let mut done = 0;
        let mut stalled = 0;
        let release_after = 10 * self.graph.qubits().max(10);
        loop {
            // Run everything that can run.
            let mut progressed = false;
            loop {
                next.clear();
                let mut ran = false;
                for &node in &front {
                    if self.can_run(node) {
                        self.emit(node);
                        done += 1;
                        ran = true;
                        for &s in dag.successors(node) {
                            remaining[s] -= 1;
                            if remaining[s] == 0 {
                                next.push(s);
                            }
                        }
                    } else {
                        next.push(node);
                    }
                }
                std::mem::swap(&mut front, &mut next);
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
            self.extended_set(&front, &remaining);
            let (a, b) = self.best_swap(&front);
            self.swap(a, b);
            stalled += 1;
            if self.swaps.is_multiple_of(5) {
                self.decay.iter_mut().for_each(|d| *d = 1.0);
            }
        }
        Trace {
            steps: self.steps,
            initial: self.initial,
            final_layout: self.physical,
            swaps: self.swaps,
            swap_cost: self.swap_cost,
        }
    }

    fn can_run(&self, node: usize) -> bool {
        match self.dag.kinds[node] {
            Kind::TwoQubit(_) => {
                let q = &self.dag.circuit[node].qubits;
                self.graph.coupled(self.physical[q[0]], self.physical[q[1]])
            }
            _ => true,
        }
    }

    fn emit(&mut self, node: usize) {
        if self.options.merge_weight > 0.0 {
            self.track(node);
        }
        self.steps.push(Step::Node(node));
    }

    /// Add node `node`, on the physical qubits it now has, to the open blocks.
    fn track(&mut self, node: usize) {
        let qubits = &self.dag.circuit[node].qubits;
        match self.dag.kinds[node] {
            Kind::OneQubit(g) => {
                let g = g.map(|i| self.dag.ones[i as usize]);
                let p = self.physical[qubits[0]];
                if let Some(b) = self.open[p] {
                    let block = &mut self.blocks[b];
                    block.matrix = match (block.matrix, g) {
                        (Some(m), Some(g)) => Some(left_1q(&g, &m, block.pair[0] == p)),
                        _ => None,
                    };
                    block.saving = None;
                }
            }
            Kind::TwoQubit(g) => {
                let g = g.map(|i| &self.dag.twos[i as usize]);
                let (p, q) = (self.physical[qubits[0]], self.physical[qubits[1]]);
                self.track_two(p, q, g);
            }
            Kind::Other => {
                for &q in qubits {
                    self.open[self.physical[q]] = None;
                }
            }
        }
    }

    /// Add a two-qubit gate on physical qubits `(p, q)` to the open blocks:
    /// `g` is its matrix on `(p, q)` and on `(q, p)`.
    fn track_two(&mut self, p: usize, q: usize, g: Option<&[M4; 2]>) {
        match (self.open[p], self.open[q]) {
            (Some(b), Some(c)) if b == c => {
                let block = &mut self.blocks[b];
                let oriented = g.map(|g| g[usize::from(block.pair[0] != p)]);
                block.matrix = match (block.matrix, oriented) {
                    (Some(m), Some(g)) => Some(g * m),
                    _ => None,
                };
                block.saving = None;
            }
            _ => {
                self.blocks.push(OpenBlock {
                    pair: [p, q],
                    matrix: g.map(|g| g[0]),
                    saving: None,
                });
                self.open[p] = Some(self.blocks.len() - 1);
                self.open[q] = Some(self.blocks.len() - 1);
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
                Some(m) => *self
                    .savings
                    .borrow_mut()
                    .entry(bits(&m))
                    .or_insert_with(|| {
                        (cx_count(&m) + 3).saturating_sub(cx_count(&(swap_matrix() * m)))
                    }),
                None => 0,
            });
        }
        block.saving.unwrap()
    }

    fn swap(&mut self, a: usize, b: usize) {
        if self.options.merge_weight > 0.0 {
            self.swap_cost += 3 - self.saving(a, b).min(3);
            let swap = swap_matrix();
            self.track_two(a, b, Some(&[swap, swap]));
        } else {
            self.swap_cost += 3;
        }
        self.steps.push(Step::Swap(a, b));
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

    /// Set `self.extended` to the qubit pairs of the next two-qubit gates
    /// after the front layer.
    fn extended_set(&mut self, front: &[usize], remaining: &[usize]) {
        let dag = self.dag;
        self.extended.clear();
        self.stamp += 1;
        let stamp = self.stamp;
        self.queue.clear();
        self.queue.extend(front.iter().copied());
        // Track the predecessors met so far, so only gates whose blocking
        // predecessors are all in the front layer or the set count early.
        while let Some(node) = self.queue.pop_front() {
            for &s in dag.successors(node) {
                if self.pending_stamp[s] != stamp {
                    self.pending_stamp[s] = stamp;
                    self.pending[s] = remaining[s];
                }
                self.pending[s] = self.pending[s].saturating_sub(1);
                if self.pending[s] == 0 && self.seen[s] != stamp {
                    self.seen[s] = stamp;
                    if let Kind::TwoQubit(_) = dag.kinds[s] {
                        let q = &dag.circuit[s].qubits;
                        self.extended.push((q[0], q[1]));
                        if self.extended.len() >= self.options.extended_set {
                            return;
                        }
                    }
                    self.queue.push_back(s);
                }
            }
        }
        let room = self.options.extended_set - self.extended.len();
        self.extended.extend(self.ahead.iter().take(room));
    }

    fn best_swap(&mut self, front: &[usize]) -> (usize, usize) {
        let dag = self.dag;
        // The gates' pairs on physical qubits, once for all candidates.
        let physical = &self.physical;
        let on_chip = |q: &[usize]| (physical[q[0]], physical[q[1]]);
        self.gates.clear();
        self.gates
            .extend(front.iter().filter_map(|&g| match dag.kinds[g] {
                Kind::TwoQubit(_) => Some(on_chip(&dag.circuit[g].qubits)),
                _ => None,
            }));
        self.extended_on_chip.clear();
        self.extended_on_chip
            .extend(self.extended.iter().map(|&(a, b)| on_chip(&[a, b])));
        self.candidates.clear();
        for &(a, b) in &self.gates {
            for p in [a, b] {
                for &nb in self.graph.neighbours(p) {
                    let edge = (p.min(nb), p.max(nb));
                    if !self.candidates.contains(&edge) {
                        self.candidates.push(edge);
                    }
                }
            }
        }
        self.best.clear();
        let mut best_score = f64::INFINITY;
        let graph = self.graph;
        let front_now: f64 = self
            .gates
            .iter()
            .map(|&(a, b)| graph.distance(a, b) as f64)
            .sum();
        for c in 0..self.candidates.len() {
            let (p1, p2) = self.candidates[c];
            let moved = |p: usize| {
                if p == p1 {
                    p2
                } else if p == p2 {
                    p1
                } else {
                    p
                }
            };
            // The mean distance after the SWAP. Distances are whole numbers,
            // so their sum is exact whichever way it is added up, and the
            // mean is the same float as a sum of floats would give.
            let cost = |pairs: &[(usize, usize)]| -> f64 {
                if pairs.is_empty() {
                    return 0.0;
                }
                let sum: u32 = pairs
                    .iter()
                    .map(|&(a, b)| graph.distance(moved(a), moved(b)))
                    .sum();
                sum as f64 / pairs.len() as f64
            };
            let front_cost = cost(&self.gates);
            let mut score = self.decay[p1].max(self.decay[p2])
                * (front_cost + self.options.extended_weight * cost(&self.extended_on_chip));
            // Merging is a bonus only for SWAPs that bring the front closer:
            // a free SWAP that does not would be chosen again and again.
            let closer = front_cost * (self.gates.len() as f64) < front_now - 1e-9;
            if self.options.merge_weight > 0.0 && closer {
                score -= self.options.merge_weight * self.saving(p1, p2) as f64;
            }
            if score < best_score - 1e-12 {
                best_score = score;
                self.best.clear();
                self.best.push((p1, p2));
            } else if (score - best_score).abs() <= 1e-12 {
                self.best.push((p1, p2));
            }
        }
        self.best[self.rng.below(self.best.len())]
    }

    /// Route the closest front gate along a shortest path.
    fn release_valve(&mut self, front: &[usize]) {
        let dag = self.dag;
        let (a, b) = front
            .iter()
            .filter(|&&g| matches!(dag.kinds[g], Kind::TwoQubit(_)))
            .map(|&g| &dag.circuit[g])
            .map(|i| (self.physical[i.qubits[0]], self.physical[i.qubits[1]]))
            .min_by_key(|&(a, b)| self.graph.distance(a, b))
            .expect("a blocked gate");
        let path = self.graph.path(a, b);
        for w in path.windows(2).take(path.len().saturating_sub(2)) {
            self.swap(w[0], w[1]);
        }
    }
}

fn swap_matrix() -> M4 {
    matrix_2q("swap", &[]).expect("a standard gate")
}

/// `(g ⊗ I) · m` if `high`, else `(I ⊗ g) · m`: a single-qubit gate after a
/// block. Each entry is the sum of the two terms of the full product that
/// are not multiplied by a zero of `g ⊗ I`, added in the same order, so the
/// result is the full product's (the routing must not depend on how it is
/// computed), at a quarter of its multiplications.
fn left_1q(g: &M2, m: &M4, high: bool) -> M4 {
    let (g, m) = (&g.0, &m.0);
    let mut out = [[crate::linalg::ZERO; 4]; 4];
    for (r, row) in out.iter_mut().enumerate() {
        // Row r = 2i + k of g ⊗ I is g[i][0], g[i][1] at columns k, 2 + k;
        // of I ⊗ g it is g[k][0], g[k][1] at columns 2i, 2i + 1.
        let (i, k) = (r / 2, r % 2);
        let (coefficients, rows) = if high {
            (g[i], [k, 2 + k])
        } else {
            (g[k], [2 * i, 2 * i + 1])
        };
        for (c, value) in row.iter_mut().enumerate() {
            *value = coefficients[0] * m[rows[0]][c] + coefficients[1] * m[rows[1]][c];
        }
    }
    M4(out)
}

/// Gate `g` with its operands the other way round: `SWAP · g · SWAP`, which
/// only exchanges the indices 01 and 10 of rows and columns.
fn reversed(g: &M4) -> M4 {
    const OTHER_WAY: [usize; 4] = [0, 2, 1, 3];
    let mut out = g.0;
    for (i, row) in out.iter_mut().enumerate() {
        for (j, value) in row.iter_mut().enumerate() {
            *value = g.0[OTHER_WAY[i]][OTHER_WAY[j]];
        }
    }
    M4(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::linalg::kron;
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

    fn better(a: &Routed, b: &Routed) -> bool {
        (a.swap_cost, a.swaps) < (b.swap_cost, b.swaps)
    }

    #[test]
    fn trials_side_by_side_choose_as_one_after_the_other() {
        // Trial t is the single trial of seed + t·0x9E3779B9, so a loop over
        // single trials is the search as it ran before the trials were
        // spread over threads: the first of the cheapest must be chosen.
        let mut rng = TestRng::new(12);
        for graph in [CouplingGraph::line(6), CouplingGraph::grid(2, 4)] {
            let n = graph.qubits();
            for seed in 0..4 {
                let circuit = random_circuit(&mut rng, n, 40);
                let options = SabreOptions {
                    trials: 6,
                    seed,
                    ..SabreOptions::default()
                };
                let mut serial: Option<Routed> = None;
                for trial in 0..options.trials {
                    let single = SabreOptions {
                        trials: 1,
                        seed: seed.wrapping_add(trial as u64 * 0x9E37_79B9),
                        ..options
                    };
                    let routed = layout_and_route(&circuit, n, &graph, &single).unwrap();
                    if serial.as_ref().is_none_or(|b| better(&routed, b)) {
                        serial = Some(routed);
                    }
                }
                let serial = serial.unwrap();
                let parallel = layout_and_route(&circuit, n, &graph, &options).unwrap();
                assert_eq!(parallel.instructions, serial.instructions);
                assert_eq!(parallel.initial_layout, serial.initial_layout);
                assert_eq!(parallel.final_layout, serial.final_layout);
            }
        }
    }

    /// Equal entry by entry, where a zero may have either sign.
    fn same_values(a: &M4, b: &M4) -> bool {
        a.0.iter()
            .flatten()
            .zip(b.0.iter().flatten())
            .all(|(x, y)| x.re == y.re && x.im == y.im)
    }

    #[test]
    fn a_single_qubit_gate_after_a_block_is_the_full_product() {
        // Routing decisions depend on these matrices through CNOT counts at
        // tolerance edges, so the shortcut must give the very same numbers
        // as the product it replaces, not merely close ones.
        let mut rng = TestRng::new(31);
        for _ in 0..50 {
            let mut block = M4::IDENTITY;
            for _ in 0..6 {
                let angles: Vec<f64> = (0..3).map(|_| rng.range(-3.0, 3.0)).collect();
                let u3 = matrix_1q("u3", &angles).unwrap();
                block = match rng.below(3) {
                    0 => kron(&u3, &M2::IDENTITY) * block,
                    1 => kron(&M2::IDENTITY, &u3) * block,
                    _ => matrix_2q("rzz", &angles[..1]).unwrap() * block,
                };
            }
            for name in ["h", "s", "t", "sx", "x", "rz", "u3"] {
                let params: Vec<f64> = match name {
                    "rz" => vec![rng.range(-3.0, 3.0)],
                    "u3" => (0..3).map(|_| rng.range(-3.0, 3.0)).collect(),
                    _ => vec![],
                };
                let g = matrix_1q(name, &params).unwrap();
                assert!(same_values(
                    &left_1q(&g, &block, true),
                    &(kron(&g, &M2::IDENTITY) * block)
                ));
                assert!(same_values(
                    &left_1q(&g, &block, false),
                    &(kron(&M2::IDENTITY, &g) * block)
                ));
            }
        }
    }

    #[test]
    fn a_gate_the_other_way_round_is_conjugated_by_swaps() {
        let swap = swap_matrix();
        let mut rng = TestRng::new(32);
        for (name, params) in [
            ("cx", vec![]),
            ("cz", vec![]),
            ("cy", vec![]),
            ("rzz", vec![0.3]),
            ("cp", vec![rng.range(-3.0, 3.0)]),
            ("crx", vec![rng.range(-3.0, 3.0)]),
        ] {
            let g = matrix_2q(name, &params).unwrap();
            assert!(same_values(&reversed(&g), &(swap * g * swap)), "{name}");
        }
    }

    #[test]
    fn the_savings_memo_changes_nothing() {
        // Routings that share a memo of savings route exactly as routings
        // that each work everything out again.
        let mut rng = TestRng::new(33);
        let graph = CouplingGraph::grid(3, 3);
        let options = SabreOptions {
            merge_weight: 0.5,
            ..SabreOptions::default()
        };
        let identity: Vec<usize> = (0..9).collect();
        for seed in 0..4 {
            let circuit = random_circuit(&mut rng, 9, 120);
            let dag = Dag::new(&circuit, 9, &options);
            let shared = Savings::default();
            for swap_trial in 0..3 {
                let options = SabreOptions {
                    seed: seed * 10 + swap_trial,
                    ..options
                };
                let with = Router::new(&dag, &shared, &graph, &identity, &options).run();
                let alone =
                    Router::new(&dag, &Savings::default(), &graph, &identity, &options).run();
                assert!(with.steps == alone.steps);
                assert_eq!(with.swap_cost, alone.swap_cost);
            }
            assert!(!shared.borrow().is_empty());
        }
    }

    #[test]
    fn the_first_swap_trials_choose_as_a_narrower_search() {
        // The second routing is the one a search with only `first` swap
        // trials per layout would have kept; none if it is the same.
        let mut rng = TestRng::new(34);
        let graph = CouplingGraph::grid(3, 3);
        let mut differed = 0;
        for seed in 0..6 {
            let circuit = random_circuit(&mut rng, 9, 150);
            let wide = SabreOptions {
                trials: 4,
                swap_trials: 12,
                merge_weight: 0.5,
                seed,
                ..SabreOptions::default()
            };
            let narrow = SabreOptions {
                swap_trials: 3,
                ..wide
            };
            let both = layout_and_route_and_first(&circuit, 9, &graph, &wide, 3).unwrap();
            let best = layout_and_route(&circuit, 9, &graph, &wide).unwrap();
            let first = layout_and_route(&circuit, 9, &graph, &narrow).unwrap();
            assert_eq!(both[0].instructions, best.instructions);
            assert!(both[0].swap_cost <= first.swap_cost);
            match both.get(1) {
                Some(second) => {
                    differed += 1;
                    assert_eq!(second.instructions, first.instructions);
                    assert_eq!(second.initial_layout, first.initial_layout);
                    check_equivalent(&circuit, second, 9);
                }
                None => assert_eq!(best.instructions, first.instructions),
            }
            assert_eq!(
                both.len(),
                1 + usize::from(
                    best.instructions != first.instructions
                        || best.initial_layout != first.initial_layout
                )
            );
        }
        assert!(differed > 0);
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
        let graph = CouplingGraph::grid(3, 3);
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
