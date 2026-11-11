//! A fast hasher for the passes' internal maps.
//!
//! The standard library's `HashMap` hashes with SipHash-1-3, which resists
//! keys chosen to collide but costs tens of instructions per word. The keys
//! here are qubit numbers, parities and the bits of small matrices, none of
//! them chosen by an adversary, and hashing them was about 4% of compiling a
//! Quantum Volume circuit. [`FxHasher`] is the multiply-and-rotate hash of
//! Firefox and rustc (`rustc-hash`): one rotate, xor and multiply per word.
//!
//! It is also deterministic, where `RandomState` is seeded per process. No
//! pass lets a map's iteration order reach its output (the benchmark
//! outputs are the same byte for byte with either hasher), and with this one
//! that could not vary between runs even if a pass did.

use std::collections::{HashMap, HashSet};
use std::hash::{BuildHasherDefault, Hasher};

/// The multiplier of rustc's FxHasher: an odd constant with well-spread bits.
const SEED: u64 = 0x51_7c_c1_b7_27_22_0a_95;

#[derive(Clone, Copy, Default)]
pub struct FxHasher {
    hash: u64,
}

impl FxHasher {
    #[inline]
    fn add(&mut self, word: u64) {
        self.hash = (self.hash.rotate_left(5) ^ word).wrapping_mul(SEED);
    }
}

impl Hasher for FxHasher {
    #[inline]
    fn write(&mut self, bytes: &[u8]) {
        let mut chunks = bytes.chunks_exact(8);
        for chunk in &mut chunks {
            self.add(u64::from_le_bytes(chunk.try_into().unwrap()));
        }
        let rest = chunks.remainder();
        if !rest.is_empty() {
            let mut word = [0u8; 8];
            word[..rest.len()].copy_from_slice(rest);
            self.add(u64::from_le_bytes(word));
        }
    }

    #[inline]
    fn write_u8(&mut self, i: u8) {
        self.add(i as u64);
    }

    #[inline]
    fn write_u32(&mut self, i: u32) {
        self.add(i as u64);
    }

    #[inline]
    fn write_u64(&mut self, i: u64) {
        self.add(i);
    }

    #[inline]
    fn write_usize(&mut self, i: usize) {
        self.add(i as u64);
    }

    #[inline]
    fn finish(&self) -> u64 {
        self.hash
    }
}

pub type FastMap<K, V> = HashMap<K, V, BuildHasherDefault<FxHasher>>;
pub type FastSet<K> = HashSet<K, BuildHasherDefault<FxHasher>>;

#[cfg(test)]
mod tests {
    use super::*;
    use std::hash::Hash;

    fn hash<T: Hash + ?Sized>(value: &T) -> u64 {
        let mut h = FxHasher::default();
        value.hash(&mut h);
        h.finish()
    }

    #[test]
    fn equal_keys_hash_equal_and_nearby_keys_apart() {
        assert_eq!(hash(&vec![1usize, 2, 3]), hash(&[1usize, 2, 3][..]));
        assert_ne!(hash(&vec![1usize, 2]), hash(&vec![2usize, 1]));
        assert_ne!(hash(&"cx"), hash(&"cz"));
        // A byte string that is not a whole number of words still counts
        // every byte.
        assert_ne!(hash(&"swap1"), hash(&"swap2"));
        let distinct: FastSet<u64> = (0..1000usize).map(|q| hash(&q)).collect();
        assert_eq!(distinct.len(), 1000);
    }

    #[test]
    fn maps_work_as_usual() {
        let mut map: FastMap<[usize; 2], f64> = FastMap::default();
        map.insert([0, 1], 0.5);
        *map.entry([0, 1]).or_insert(0.0) += 0.25;
        map.insert([1, 0], 1.0);
        assert_eq!(map[&[0, 1]], 0.75);
        assert_eq!(map.len(), 2);
    }
}
