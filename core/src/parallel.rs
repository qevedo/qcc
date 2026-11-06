//! Independent work side by side, with results that do not depend on the
//! threads.
//!
//! Much of compilation is several candidates computed from the same input
//! and compared at the end: SABRE's layout trials, the direct pipeline and
//! the one through CX for CZ devices, the logical circuit with and without a
//! Pauli network on a chip. None reads another's result, so they can run on
//! different cores. Each helper here returns the results in the order of the
//! work, as a loop would, so that the caller's choice among them (the first
//! of the cheapest, the first error) is the same as before and the output
//! does not depend on which thread finished first. On WebAssembly, which has
//! no threads, the work runs one after the other.

use std::sync::atomic::{AtomicUsize, Ordering};

/// `f(0)`, …, `f(n - 1)`, in that order, computed on up to as many threads
/// as there are cores. Each thread takes the next index until none is left,
/// as the costs of the items differ.
pub fn map_ordered<R: Send>(n: usize, f: impl Fn(usize) -> R + Sync) -> Vec<R> {
    let threads = std::thread::available_parallelism().map_or(1, |t| t.get());
    if n < 2 || threads < 2 || cfg!(target_arch = "wasm32") {
        return (0..n).map(f).collect();
    }
    let next = AtomicUsize::new(0);
    let (f, next) = (&f, &next);
    let mut done: Vec<(usize, R)> = std::thread::scope(|scope| {
        let handles: Vec<_> = (0..threads.min(n))
            .map(|_| {
                scope.spawn(move || {
                    let mut out = Vec::new();
                    loop {
                        let i = next.fetch_add(1, Ordering::Relaxed);
                        if i >= n {
                            break;
                        }
                        out.push((i, f(i)));
                    }
                    out
                })
            })
            .collect();
        handles
            .into_iter()
            .flat_map(|h| h.join().expect("a worker thread"))
            .collect()
    });
    done.sort_unstable_by_key(|&(i, _)| i);
    done.into_iter().map(|(_, r)| r).collect()
}

/// [`join`] if `parallel`, else `a()` then `b()` on this thread: for work
/// too small to pay for a thread.
pub fn join_if<A: Send, B: Send>(
    parallel: bool,
    a: impl FnOnce() -> A + Send,
    b: impl FnOnce() -> B + Send,
) -> (A, B) {
    if parallel { join(a, b) } else { (a(), b()) }
}

/// `a()` and `b()`, `b` on a second thread while `a` runs on this one.
pub fn join<A: Send, B: Send>(
    a: impl FnOnce() -> A + Send,
    b: impl FnOnce() -> B + Send,
) -> (A, B) {
    if cfg!(target_arch = "wasm32") {
        return (a(), b());
    }
    std::thread::scope(|scope| {
        let b = scope.spawn(b);
        let a = a();
        (a, b.join().expect("a worker thread"))
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn results_come_in_the_order_of_the_work() {
        // Uneven costs, so that threads finish out of order.
        let out = map_ordered(64, |i| {
            let mut x = i as u64;
            for _ in 0..(64 - i) * 2000 {
                x = x.wrapping_mul(6364136223846793005).wrapping_add(1);
            }
            (i, x)
        });
        assert_eq!(out.len(), 64);
        for (k, &(i, _)) in out.iter().enumerate() {
            assert_eq!(i, k);
        }
        assert!(map_ordered(0, |i| i).is_empty());
        assert_eq!(map_ordered(1, |i| i + 7), vec![7]);
    }

    #[test]
    fn join_returns_both_results_in_place() {
        let (a, b) = join(|| "a".repeat(3), || (1..=10).sum::<i32>());
        assert_eq!((a.as_str(), b), ("aaa", 55));
    }
}
