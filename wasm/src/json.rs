//! Writing the JSON replies by hand: the core has no JSON dependency.

pub fn error(message: &str) -> String {
    format!(r#"{{"ok":false,"message":{}}}"#, string(message))
}

pub fn number(x: f64) -> String {
    if x.is_finite() {
        format!("{x:?}")
    } else {
        "null".into()
    }
}

pub fn string(s: &str) -> String {
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
