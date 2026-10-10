//! Lowering every gate to a device's native gates with the fewest gates.
//! Matches `qevedo/compiler/passes/decompose.py`.

use crate::gates::{is_known_gate, matrix_1q, matrix_2q};
use crate::ir::{Instruction, NON_GATES, cost};
use crate::linalg::{M2, M4};
use crate::one_qubit::{Basis1q, bases_for, count_1q, synthesize_1q};
use crate::passes::merge_single_qubit_runs;
use crate::two_qubit::{Step, TwoQubitBasis, TwoQubitSynth, bases_2q};
use std::collections::{HashMap, HashSet};
use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI};

/// A gate of a definition: name, parameters and operand indices, in time order.
type Body = Vec<(&'static str, Vec<f64>, Vec<usize>)>;

fn g(
    name: &'static str,
    params: &[f64],
    operands: &[usize],
) -> (&'static str, Vec<f64>, Vec<usize>) {
    (name, params.to_vec(), operands.to_vec())
}

fn crz(t: f64) -> Body {
    vec![
        g("rz", &[t / 2.0], &[1]),
        g("cx", &[], &[0, 1]),
        g("rz", &[-t / 2.0], &[1]),
        g("cx", &[], &[0, 1]),
    ]
}

fn cp(t: f64) -> Body {
    vec![
        g("p", &[t / 2.0], &[0]),
        g("cx", &[], &[0, 1]),
        g("p", &[-t / 2.0], &[1]),
        g("cx", &[], &[0, 1]),
        g("p", &[t / 2.0], &[1]),
    ]
}

fn cu(theta: f64, phi: f64, lam: f64, gamma: f64) -> Body {
    vec![
        g("p", &[gamma + (lam + phi) / 2.0], &[0]),
        g("p", &[(lam - phi) / 2.0], &[1]),
        g("cx", &[], &[0, 1]),
        g("u3", &[-theta / 2.0, 0.0, -(phi + lam) / 2.0], &[1]),
        g("cx", &[], &[0, 1]),
        g("u3", &[theta / 2.0, phi, 0.0], &[1]),
    ]
}

fn rzz(t: f64) -> Body {
    vec![
        g("cx", &[], &[0, 1]),
        g("rz", &[t], &[1]),
        g("cx", &[], &[0, 1]),
    ]
}

fn wrapped(before: Body, body: Body, after: Body) -> Body {
    before.into_iter().chain(body).chain(after).collect()
}

/// The textbook definition of a standard gate in smaller gates.
pub fn definition(name: &str, p: &[f64]) -> Option<Body> {
    let h1 = || vec![g("h", &[], &[1])];
    Some(match (name, p) {
        ("swap", []) => vec![
            g("cx", &[], &[0, 1]),
            g("cx", &[], &[1, 0]),
            g("cx", &[], &[0, 1]),
        ],
        ("cz", []) => wrapped(h1(), vec![g("cx", &[], &[0, 1])], h1()),
        ("cy", []) => vec![
            g("sdg", &[], &[1]),
            g("cx", &[], &[0, 1]),
            g("s", &[], &[1]),
        ],
        ("ch", []) => vec![
            g("ry", &[FRAC_PI_4], &[1]),
            g("cx", &[], &[0, 1]),
            g("ry", &[-FRAC_PI_4], &[1]),
        ],
        ("csx", []) => wrapped(h1(), cp(FRAC_PI_2), h1()),
        ("crz", [t]) => crz(*t),
        ("crx", [t]) => wrapped(h1(), crz(*t), h1()),
        ("cry", [t]) => vec![
            g("ry", &[t / 2.0], &[1]),
            g("cx", &[], &[0, 1]),
            g("ry", &[-t / 2.0], &[1]),
            g("cx", &[], &[0, 1]),
        ],
        ("cp" | "cphase" | "cu1", [t]) => cp(*t),
        ("cu3", [a, b, c]) => cu(*a, *b, *c, 0.0),
        ("cu", [a, b, c, d]) => cu(*a, *b, *c, *d),
        ("rzz", [t]) => rzz(*t),
        ("rxx", [t]) => wrapped(
            vec![g("h", &[], &[0]), g("h", &[], &[1])],
            rzz(*t),
            vec![g("h", &[], &[0]), g("h", &[], &[1])],
        ),
        ("ryy", [t]) => wrapped(
            vec![g("rx", &[FRAC_PI_2], &[0]), g("rx", &[FRAC_PI_2], &[1])],
            rzz(*t),
            vec![g("rx", &[-FRAC_PI_2], &[0]), g("rx", &[-FRAC_PI_2], &[1])],
        ),
        ("rzx", [t]) => wrapped(h1(), rzz(*t), h1()),
        ("dcx", []) => vec![g("cx", &[], &[0, 1]), g("cx", &[], &[1, 0])],
        ("iswap", []) => vec![
            g("s", &[], &[0]),
            g("s", &[], &[1]),
            g("h", &[], &[0]),
            g("cx", &[], &[0, 1]),
            g("cx", &[], &[1, 0]),
            g("h", &[], &[1]),
        ],
        // Toffoli with 6 CX, which is optimal.
        ("ccx", []) => vec![
            g("h", &[], &[2]),
            g("cx", &[], &[1, 2]),
            g("tdg", &[], &[2]),
            g("cx", &[], &[0, 2]),
            g("t", &[], &[2]),
            g("cx", &[], &[1, 2]),
            g("tdg", &[], &[2]),
            g("cx", &[], &[0, 2]),
            g("t", &[], &[1]),
            g("t", &[], &[2]),
            g("h", &[], &[2]),
            g("cx", &[], &[0, 1]),
            g("t", &[], &[0]),
            g("tdg", &[], &[1]),
            g("cx", &[], &[0, 1]),
        ],
        ("cswap", []) => vec![
            g("cx", &[], &[2, 1]),
            g("ccx", &[], &[0, 1, 2]),
            g("cx", &[], &[2, 1]),
        ],
        // Toffoli up to a relative phase, with 3 CX.
        ("rccx", []) => vec![
            g("u2", &[0.0, PI], &[2]),
            g("u1", &[FRAC_PI_4], &[2]),
            g("cx", &[], &[1, 2]),
            g("u1", &[-FRAC_PI_4], &[2]),
            g("cx", &[], &[0, 2]),
            g("u1", &[FRAC_PI_4], &[2]),
            g("cx", &[], &[1, 2]),
            g("u1", &[-FRAC_PI_4], &[2]),
            g("u2", &[0.0, PI], &[2]),
        ],
        _ => return None,
    })
}

/// Lowers instructions to a fixed native gate set, remembering every gate it lowered.
pub struct Lowering {
    pub native: HashSet<String>,
    pub bases_1q: Vec<Basis1q>,
    bases_2q: Vec<TwoQubitBasis>,
    synth: TwoQubitSynth,
    memo: HashMap<(String, Vec<u64>, usize), Vec<Instruction>>,
    memo_2q: HashMap<[u64; 32], Vec<Instruction>>,
}

impl Lowering {
    pub fn new(native: impl IntoIterator<Item = String>) -> Lowering {
        let native: HashSet<String> = native.into_iter().map(|g| g.to_lowercase()).collect();
        Lowering {
            bases_1q: bases_for(&native),
            bases_2q: bases_2q(&native),
            native,
            synth: TwoQubitSynth::new(),
            memo: HashMap::new(),
            memo_2q: HashMap::new(),
        }
    }

    /// `inst` in native gates.
    pub fn lower(&mut self, inst: &Instruction) -> Result<Vec<Instruction>, String> {
        let name = inst.name.to_lowercase();
        if NON_GATES.contains(&name.as_str()) || self.native.contains(&name) {
            return Ok(vec![inst.clone()]);
        }
        // Lower each distinct gate once, on placeholder qubits, then relabel.
        let key = (
            name.clone(),
            inst.params.iter().map(|p| p.to_bits()).collect(),
            inst.qubits.len(),
        );
        if !self.memo.contains_key(&key) {
            let placeholder =
                Instruction::gate(&name, (0..inst.qubits.len()).collect(), inst.params.clone());
            let lowered = self.lower_new(&placeholder)?;
            self.memo.insert(key.clone(), lowered);
        }
        Ok(self.memo[&key]
            .iter()
            .map(|low| low.relabel(&inst.qubits))
            .collect())
    }

    fn lower_new(&mut self, inst: &Instruction) -> Result<Vec<Instruction>, String> {
        let name = inst.name.as_str();
        let mut candidates: Vec<Vec<Instruction>> = Vec::new();
        if let Some(body) = definition(name, &inst.params) {
            let mut out = Vec::new();
            for (sub, params, operands) in body {
                let qubits = operands.iter().map(|&i| inst.qubits[i]).collect();
                out.extend(self.lower(&Instruction::gate(sub, qubits, params))?);
            }
            candidates.push(merge_single_qubit_runs(out, &self.native, &self.bases_1q));
        }
        if is_known_gate(name) && inst.qubits.len() <= 2 {
            if inst.qubits.len() == 1 {
                let m = matrix_1q(name, &inst.params).ok_or_else(|| bad_params(inst))?;
                candidates.push(self.one_qubit(&m, inst.qubits[0])?);
            } else if !self.bases_2q.is_empty() {
                let m = matrix_2q(name, &inst.params).ok_or_else(|| bad_params(inst))?;
                candidates.push(self.two_qubit(&m, [inst.qubits[0], inst.qubits[1]])?);
            }
        }
        if candidates.is_empty() {
            if !is_known_gate(name) && definition(name, &inst.params).is_none() {
                return Err(format!(
                    "cannot lower gate '{name}': qcc does not know its definition yet"
                ));
            }
            let mut native: Vec<&String> = self.native.iter().collect();
            native.sort();
            return Err(format!(
                "cannot lower '{name}': the native gates {native:?} include no supported two-qubit gate \
                 (cx, cz, cy, ch, ecr, rxx, ryy, rzz or rzx)"
            ));
        }
        // The first of the cheapest, like Python's min.
        let mut best = candidates.swap_remove(0);
        for candidate in candidates {
            if cost(&candidate) < cost(&best) {
                best = candidate;
            }
        }
        Ok(best)
    }

    fn one_qubit(&self, m: &M2, qubit: usize) -> Result<Vec<Instruction>, String> {
        if self.bases_1q.is_empty() {
            let mut native: Vec<&String> = self.native.iter().collect();
            native.sort();
            return Err(format!(
                "the native gates {native:?} have no universal single-qubit basis"
            ));
        }
        Ok(synthesize_1q(m, &self.bases_1q)
            .into_iter()
            .map(|op| Instruction::gate(&op.name, vec![qubit], op.params))
            .collect())
    }

    /// The cheapest native sequence for the two-qubit unitary `m` on `qubits`.
    pub fn two_qubit(&mut self, m: &M4, qubits: [usize; 2]) -> Result<Vec<Instruction>, String> {
        let mut key = [0u64; 32];
        for i in 0..4 {
            for j in 0..4 {
                key[8 * i + 2 * j] = m.0[i][j].re.to_bits();
                key[8 * i + 2 * j + 1] = m.0[i][j].im.to_bits();
            }
        }
        if !self.memo_2q.contains_key(&key) {
            let lowered = self.two_qubit_new(m)?;
            self.memo_2q.insert(key, lowered);
        }
        Ok(self.memo_2q[&key]
            .iter()
            .map(|low| low.relabel(&qubits))
            .collect())
    }

    fn two_qubit_new(&mut self, m: &M4) -> Result<Vec<Instruction>, String> {
        if self.bases_1q.is_empty() {
            self.one_qubit(&M2::IDENTITY, 0)?;
        }
        let bases_1q = &self.bases_1q;
        let cost_1q = |u: &M2| {
            if bases_1q.is_empty() {
                0
            } else {
                count_1q(u, bases_1q)
            }
        };
        let mut best: Option<Vec<Instruction>> = None;
        for basis in &self.bases_2q {
            let mut out = Vec::new();
            for step in self.synth.synthesize(m, basis, Some(&cost_1q)) {
                match step {
                    Step::U(q, u) => out.extend(
                        synthesize_1q(&u, bases_1q)
                            .into_iter()
                            .map(|op| Instruction::gate(&op.name, vec![q], op.params)),
                    ),
                    Step::Gate(name, params, q) => {
                        out.push(Instruction::gate(&name, q.to_vec(), params))
                    }
                }
            }
            if best.as_ref().is_none_or(|b| cost(&out) < cost(b)) {
                best = Some(out);
            }
        }
        Ok(best.expect("two_qubit needs a two-qubit basis"))
    }
}

fn bad_params(inst: &Instruction) -> String {
    format!(
        "gate '{}' does not take {} parameter(s)",
        inst.name,
        inst.params.len()
    )
}
