//! The QCC core compiled to WebAssembly for the playgrounds on qevedo.com.
//!
//! No bindings generator: the module exports `alloc`, `dealloc` and one
//! function per tool, `decompose` and `map`. JavaScript writes a request into
//! memory from `alloc`, calls the tool with `(ptr, len)`, and reads the reply
//! at the returned pointer: a little-endian u32 length followed by that many
//! bytes of JSON. It then frees the reply with `dealloc(ptr, length + 4)`.
//!
//! | file | tool |
//! |---|---|
//! | `decompose.rs` | the gate decomposition playground |
//! | `multicore.rs` | the multi-core mapping playground |
//! | `json.rs` | writing the JSON replies |

mod decompose;
mod json;
mod multicore;

pub use decompose::respond as decompose_request;
pub use multicore::respond as map_request;

#[unsafe(no_mangle)]
pub extern "C" fn alloc(len: usize) -> *mut u8 {
    let mut buffer = Vec::<u8>::with_capacity(len);
    let ptr = buffer.as_mut_ptr();
    std::mem::forget(buffer);
    ptr
}

/// # Safety
/// `ptr` and `len` must come from `alloc(len)` or a tool's reply.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn dealloc(ptr: *mut u8, len: usize) {
    drop(unsafe { Vec::from_raw_parts(ptr, 0, len) });
}

/// Decompose a gate or a unitary (see `decompose.rs`).
///
/// # Safety
/// `ptr` must point to `len` readable bytes.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn decompose(ptr: *const u8, len: usize) -> *mut u8 {
    unsafe { serve(ptr, len, decompose::respond) }
}

/// Map a circuit onto a multi-core machine (see `multicore.rs`).
///
/// # Safety
/// `ptr` must point to `len` readable bytes.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn map(ptr: *const u8, len: usize) -> *mut u8 {
    unsafe { serve(ptr, len, multicore::respond) }
}

/// Read the request at `ptr`, answer it with `respond`, and return the reply
/// as a u32 length followed by the JSON.
unsafe fn serve(ptr: *const u8, len: usize, respond: fn(&str) -> String) -> *mut u8 {
    let request = unsafe { std::slice::from_raw_parts(ptr, len) };
    let reply = match std::str::from_utf8(request) {
        Ok(text) => respond(text),
        Err(_) => json::error("the request is not UTF-8"),
    };
    let mut out = Vec::with_capacity(reply.len() + 4);
    out.extend_from_slice(&(reply.len() as u32).to_le_bytes());
    out.extend_from_slice(reply.as_bytes());
    let ptr = out.as_mut_ptr();
    std::mem::forget(out);
    ptr
}
