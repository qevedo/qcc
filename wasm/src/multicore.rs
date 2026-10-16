//! The multi-core mapping playground: map a circuit's two-qubit gates onto
//! cores linked by quantum links.
//!
//! The request is lines of text:
//!
//! ```text
//! 40 4 10               circuit qubits, cores, qubits per core
//! line                  links: all-to-all, line, ring or grid:RxC
//! hqa                   method: hqa, hqa-latency, fgp-roee or naive
//! 1                     channels per link (EPR pairs at once), for the duration
//! 0 1,2 3,1 2           the two-qubit gates, in program order
//! ```
//!
//! The reply has the number of teleportations, the duration under the
//! default timing (`qcc::multicore::Timing`), the links, and where every
//! qubit is before the circuit and during each timeslice.

use crate::json::{error, number};
use qcc::multicore::{Architecture, Method, Timeslices, Timing, map_timed, mapping_duration};

/// Handle one request, as JSON.
pub fn respond(request: &str) -> String {
    match answer(request) {
        Ok(reply) => reply,
        Err(message) => error(&message),
    }
}

fn answer(request: &str) -> Result<String, String> {
    let lines: Vec<&str> = request.lines().collect();
    if lines.len() < 5 {
        return Err("incomplete request".into());
    }
    let sizes: Vec<usize> = lines[0]
        .split_whitespace()
        .map(str::parse)
        .collect::<Result<_, _>>()
        .map_err(|_| "the sizes are not whole numbers")?;
    let &[qubits, cores, per_core] = sizes.as_slice() else {
        return Err("give qubits, cores and qubits per core".into());
    };
    let arch = architecture(lines[1].trim(), cores, per_core)?;
    let method = Method::parse(lines[2].trim())?;
    let timing = Timing {
        link_capacity: lines[3].trim().parse().map_err(|_| "bad channel count")?,
        ..Timing::default()
    };
    let pairs = lines[4]
        .split(',')
        .filter(|g| !g.trim().is_empty())
        .map(|g| {
            let qs: Vec<usize> = g
                .split_whitespace()
                .map(str::parse)
                .collect::<Result<_, _>>()
                .map_err(|_| format!("bad gate '{g}'"))?;
            match qs.as_slice() {
                &[a, b] => Ok((a, b)),
                _ => Err(format!("bad gate '{g}'")),
            }
        })
        .collect::<Result<Vec<(usize, usize)>, String>>()?;
    let slices = Timeslices::from_pairs(qubits, &pairs)?.fitted_to(&arch);
    let mapping = map_timed(&slices, &arch, method, &timing)?;
    let list = |v: &[usize]| v.iter().map(usize::to_string).collect::<Vec<_>>().join(",");
    let assignments: Vec<String> = mapping
        .slices
        .iter()
        .map(|a| format!("[{}]", list(a)))
        .collect();
    let gates: Vec<String> = (0..slices.len())
        .map(|t| {
            let g: Vec<String> = slices
                .slice(t)
                .iter()
                .map(|(a, b)| format!("[{a},{b}]"))
                .collect();
            format!("[{}]", g.join(","))
        })
        .collect();
    let links: Vec<String> = arch
        .links()
        .iter()
        .map(|(a, b)| format!("[{a},{b}]"))
        .collect();
    Ok(format!(
        r#"{{"ok":true,"communications":{},"moved":{},"duration":{},"links":[{}],"initial":[{}],"slices":[{}],"gates":[{}]}}"#,
        mapping.communications(),
        mapping.moved_qubits(),
        number(mapping_duration(&mapping, &slices, &timing)),
        links.join(","),
        list(&mapping.initial),
        assignments.join(","),
        gates.join(",")
    ))
}

fn architecture(links: &str, cores: usize, per_core: usize) -> Result<Architecture, String> {
    match links {
        "all-to-all" => Architecture::new(cores, per_core),
        "line" => Architecture::line(cores, per_core),
        "ring" => Architecture::ring(cores, per_core),
        grid => {
            let shape = grid
                .strip_prefix("grid:")
                .ok_or_else(|| format!("unknown links '{grid}'"))?;
            let (rows, columns) = shape
                .split_once('x')
                .and_then(|(r, c)| Some((r.parse().ok()?, c.parse().ok()?)))
                .ok_or_else(|| format!("bad grid '{shape}'"))?;
            if rows * columns != cores {
                return Err(format!(
                    "a {rows}x{columns} grid has {} cores",
                    rows * columns
                ));
            }
            Architecture::grid(rows, columns, per_core)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::respond;

    #[test]
    fn replies() {
        let ok = respond("4 2 2\nall-to-all\nhqa\n1\n0 2,1 3,0 1\n");
        assert!(ok.starts_with(r#"{"ok":true,"communications":"#), "{ok}");
        assert!(ok.contains(r#""links":[[0,1]]"#), "{ok}");
        let grid = respond("8 4 2\ngrid:2x2\nfgp-roee\n2\n0 7,1 6,2 5\n");
        assert!(grid.contains(r#""ok":true"#), "{grid}");
        for (bad, message) in [
            ("9 2 2\nline\nhqa\n1\n0 1\n", "do not fit"),
            ("4 2 2\nstar\nhqa\n1\n0 1\n", "unknown links"),
            ("4 2 2\nline\nmagic\n1\n0 1\n", "unknown mapping method"),
            ("4 2 2\nline\nhqa\n1\n0 x\n", "bad gate"),
        ] {
            let reply = respond(bad);
            assert!(
                reply.contains(r#""ok":false"#) && reply.contains(message),
                "{reply}"
            );
        }
    }
}
