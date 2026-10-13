//! The QCC core compiled to WebAssembly for the decomposition playground.
//!
//! No bindings generator: the module exports `alloc`, `dealloc` and
//! `decompose`. JavaScript writes a request into memory from `alloc`, calls
//! `decompose(ptr, len)`, and reads the reply at the returned pointer: a
//! little-endian u32 length followed by that many bytes of JSON. It then frees
//! the reply with `dealloc(ptr, length + 4)`.
//!
//! The request is lines of text:
//!
//! ```text
//! rz,sx,x,cx            native gates
//! gate                  or: unitary
//! cu                    gate name            | 2 or 4 (the dimension)
//! 0.3,0.1,0.2,0         parameters (maybe empty) | re,im,re,im,… row by row
//! ```

use qcc::decompose::{Decomposition, decompose_gate, decompose_unitary};
use qcc::linalg::C;

#[unsafe(no_mangle)]
pub extern "C" fn alloc(len: usize) -> *mut u8 {
    let mut buffer = Vec::<u8>::with_capacity(len);
    let ptr = buffer.as_mut_ptr();
    std::mem::forget(buffer);
    ptr
}

/// # Safety
/// `ptr` and `len` must come from `alloc(len)` or a `decompose` reply.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn dealloc(ptr: *mut u8, len: usize) {
    drop(unsafe { Vec::from_raw_parts(ptr, 0, len) });
}

/// # Safety
/// `ptr` must point to `len` readable bytes.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn decompose(ptr: *const u8, len: usize) -> *mut u8 {
    let request = unsafe { std::slice::from_raw_parts(ptr, len) };
    let reply = match std::str::from_utf8(request) {
        Ok(text) => respond(text),
        Err(_) => error("the request is not UTF-8"),
    };
    let mut out = Vec::with_capacity(reply.len() + 4);
    out.extend_from_slice(&(reply.len() as u32).to_le_bytes());
    out.extend_from_slice(reply.as_bytes());
    let ptr = out.as_mut_ptr();
    std::mem::forget(out);
    ptr
}

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

fn error(message: &str) -> String {
    format!(r#"{{"ok":false,"message":{}}}"#, string(message))
}

fn number(x: f64) -> String {
    if x.is_finite() {
        format!("{x:?}")
    } else {
        "null".into()
    }
}

fn string(s: &str) -> String {
    let mut out = String::with_capacity(s.len() + 2);
    out.push('"');
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            c if (c as u32) < 0x20 => out.push_str(&format!("\\u{:04x}", c as u32)),
            c => out.push(c),
        }
    }
    out.push('"');
    out
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
