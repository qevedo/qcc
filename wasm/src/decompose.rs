//! The gate decomposition playground: decompose a gate or a unitary into a
//! device's native gates.
//!
//! The request is lines of text:
//!
//! ```text
//! rz,sx,x,cx            native gates
//! gate                  or: unitary
//! cu                    gate name            | 2 or 4 (the dimension)
//! 0.3,0.1,0.2,0         parameters (maybe empty) | re,im,re,im,… row by row
//! ```

use crate::json::{error, number, string};
use qcc::decompose::{Decomposition, decompose_gate, decompose_unitary};
use qcc::linalg::C;

/// Handle one request, as JSON.
pub fn respond(request: &str) -> String {
    let lines: Vec<&str> = request.lines().collect();
    if lines.len() < 4 {
        return error("incomplete request");
    }
    let native: Vec<String> = lines[0]
        .split(',')
        .map(|g| g.trim().to_lowercase())
        .filter(|g| !g.is_empty())
        .collect();
    let numbers = match lines[3]
        .split(',')
        .map(str::trim)
        .filter(|v| !v.is_empty())
        .map(str::parse::<f64>)
        .collect::<Result<Vec<f64>, _>>()
    {
        Ok(values) => values,
        Err(_) => return error("a number could not be read"),
    };
    let result = match lines[1].trim() {
        "gate" => decompose_gate(lines[2].trim(), &numbers, &native),
        "unitary" => {
            let dim: usize = lines[2].trim().parse().unwrap_or(0);
            if numbers.len() != 2 * dim * dim {
                return error(&format!(
                    "a {dim}×{dim} matrix needs {} numbers",
                    2 * dim * dim
                ));
            }
            let matrix = (0..dim)
                .map(|i| {
                    (0..dim)
                        .map(|j| C::new(numbers[2 * (i * dim + j)], numbers[2 * (i * dim + j) + 1]))
                        .collect()
                })
                .collect();
            decompose_unitary(&matrix, &native)
        }
        other => return error(&format!("unknown request kind '{other}'")),
    };
    match result {
        Ok(d) => success(&d),
        Err(message) => error(&message),
    }
}

fn success(d: &Decomposition) -> String {
    let gates: Vec<String> = d
        .instructions
        .iter()
        .map(|i| {
            let qubits: Vec<String> = i.qubits.iter().map(usize::to_string).collect();
            let params: Vec<String> = i.params.iter().map(|p| number(*p)).collect();
            format!(
                r#"{{"name":{},"qubits":[{}],"params":[{}]}}"#,
                string(&i.name),
                qubits.join(","),
                params.join(",")
            )
        })
        .collect();
    format!(
        r#"{{"ok":true,"qubits":{},"error":{},"gates":[{}]}}"#,
        d.qubits,
        number(d.error),
        gates.join(",")
    )
}

#[cfg(test)]
mod tests {
    use super::respond;

    #[test]
    fn replies() {
        let swap = respond("rz,sx,x,cx\ngate\nswap\n\n");
        assert!(swap.starts_with(r#"{"ok":true,"qubits":2"#), "{swap}");
        assert_eq!(swap.matches(r#""name":"cx""#).count(), 3);
        let h = respond(
            "rz,sx,x,cx\nunitary\n2\n0.7071067811865476,0,0.7071067811865476,0,0.7071067811865476,0,-0.7071067811865476,0\n",
        );
        assert!(h.contains(r#""name":"sx""#), "{h}");
        let bad = respond("rz,sx,cx\ngate\nrx\n\n");
        assert!(
            bad.contains(r#""ok":false"#) && bad.contains("takes 1 parameter"),
            "{bad}"
        );
    }
}
