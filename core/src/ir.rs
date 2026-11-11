//! The instruction list the passes work on: gates, measurements, resets and
//! barriers on numbered qubits and classical bits.
//!
//! The passes copy, relabel and rebuild instructions constantly: compiling
//! the 32-qubit Quantum Volume circuit copies or relabels about 280,000 of
//! them. With a `String` name and `Vec` qubits and parameters, each of those
//! cost two or three heap allocations, and allocating and freeing were 27%
//! of the core's instructions (10% since). So an instruction owns no heap
//! memory in the usual case: its name is a `&'static str` for every gate
//! name the compiler knows ([`Name`]), and its qubits and parameters are
//! stored inline up to three ([`Inline`]), spilling to a `Vec` only beyond
//! that (a wide barrier, or a gate the input defines). The outputs are the
//! same byte for byte as with the old representation.

use std::borrow::{Borrow, Cow};
use std::fmt;
use std::hash::{Hash, Hasher};
use std::ops::{Deref, DerefMut};

#[derive(Clone, Debug, PartialEq)]
pub struct Instruction {
    pub name: Name,
    pub qubits: Inline<usize, 3>,
    pub params: Inline<f64, 3>,
    pub clbits: Vec<usize>,
}

/// Operations that are not gates: they pass through unchanged and block optimizations.
pub const NON_GATES: [&str; 3] = ["barrier", "measure", "reset"];

/// Every name the compiler itself reads or writes. A name in this list costs
/// no allocation; any other (an opaque gate passed through) is kept as a
/// `String`. The lookup is a `match`, which compiles to a switch on the
/// length and a few comparisons.
macro_rules! known_names {
    ($($name:literal),* $(,)?) => {
        #[cfg(test)]
        const KNOWN_NAMES: &[&str] = &[$($name),*];

        fn known_name(name: &str) -> Option<&'static str> {
            match name {
                $($name => Some($name),)*
                _ => None,
            }
        }
    };
}

known_names!(
    "barrier", "measure", "reset", "id", "x", "y", "z", "h", "s", "sdg", "t", "tdg", "sx", "sxdg",
    "rx", "ry", "rz", "p", "u", "u1", "u2", "u3", "r", "cx", "cy", "cz", "ch", "swap", "iswap",
    "ecr", "rxx", "ryy", "rzz", "rzx", "crx", "cry", "crz", "cp", "cu1", "cu3", "csx", "ccx",
    "ccz", "cswap", "cnot", "cphase", "phase", "cu", "dcx", "rccx",
);

/// An instruction's name: borrowed from the known names when it is one of
/// them, so that copying it is free, and owned otherwise. It compares,
/// hashes and prints as the string it holds.
#[derive(Clone, Eq)]
pub struct Name(Cow<'static, str>);

impl Name {
    pub fn as_str(&self) -> &str {
        &self.0
    }
}

impl From<&str> for Name {
    fn from(name: &str) -> Name {
        match known_name(name) {
            Some(known) => Name(Cow::Borrowed(known)),
            None => Name(Cow::Owned(name.to_string())),
        }
    }
}

impl From<String> for Name {
    fn from(name: String) -> Name {
        match known_name(&name) {
            Some(known) => Name(Cow::Borrowed(known)),
            None => Name(Cow::Owned(name)),
        }
    }
}

impl From<&String> for Name {
    fn from(name: &String) -> Name {
        Name::from(name.as_str())
    }
}

impl From<Name> for String {
    fn from(name: Name) -> String {
        name.0.into_owned()
    }
}

impl Deref for Name {
    type Target = str;
    fn deref(&self) -> &str {
        &self.0
    }
}

impl AsRef<str> for Name {
    fn as_ref(&self) -> &str {
        &self.0
    }
}

impl Borrow<str> for Name {
    fn borrow(&self) -> &str {
        &self.0
    }
}

impl PartialEq for Name {
    fn eq(&self, other: &Name) -> bool {
        *self.0 == *other.0
    }
}

impl PartialEq<str> for Name {
    fn eq(&self, other: &str) -> bool {
        *self.0 == *other
    }
}

impl PartialEq<&str> for Name {
    fn eq(&self, other: &&str) -> bool {
        *self.0 == **other
    }
}

impl PartialEq<String> for Name {
    fn eq(&self, other: &String) -> bool {
        *self.0 == **other
    }
}

impl Hash for Name {
    fn hash<H: Hasher>(&self, state: &mut H) {
        (*self.0).hash(state)
    }
}

impl fmt::Debug for Name {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Debug::fmt(&*self.0, f)
    }
}

impl fmt::Display for Name {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

/// A short list stored inline, up to `N` items, and in a `Vec` beyond.
/// It reads as a slice (`Deref<Target = [T]>`), and compares, hashes and
/// prints as one, so it can stand where a `Vec` stood. Only for `Copy`
/// items, which is all an instruction holds (qubit numbers and angles).
#[derive(Clone)]
pub enum Inline<T: Copy + Default, const N: usize> {
    Inline { len: u8, items: [T; N] },
    Heap(Vec<T>),
}

impl<T: Copy + Default, const N: usize> Inline<T, N> {
    pub const fn new() -> Self {
        // `T::default()` is not const, so the inline items cannot be built
        // here. An empty `Vec` does not allocate either, and the first
        // `push` moves the list inline.
        Inline::Heap(Vec::new())
    }

    pub fn push(&mut self, item: T) {
        match self {
            Inline::Inline { len, items } if (*len as usize) < N => {
                items[*len as usize] = item;
                *len += 1;
            }
            Inline::Inline { len, items } => {
                let mut spilled = Vec::with_capacity(2 * N);
                spilled.extend_from_slice(&items[..*len as usize]);
                spilled.push(item);
                *self = Inline::Heap(spilled);
            }
            Inline::Heap(v) if v.is_empty() && v.capacity() == 0 => {
                let mut items = [T::default(); N];
                items[0] = item;
                *self = Inline::Inline { len: 1, items };
            }
            Inline::Heap(v) => v.push(item),
        }
    }

    pub fn as_slice(&self) -> &[T] {
        self
    }

    pub fn from_slice(slice: &[T]) -> Self {
        if slice.len() <= N {
            let mut items = [T::default(); N];
            items[..slice.len()].copy_from_slice(slice);
            Inline::Inline {
                len: slice.len() as u8,
                items,
            }
        } else {
            Inline::Heap(slice.to_vec())
        }
    }
}

impl<T: Copy + Default, const N: usize> Default for Inline<T, N> {
    fn default() -> Self {
        Inline::new()
    }
}

impl<T: Copy + Default, const N: usize> Deref for Inline<T, N> {
    type Target = [T];
    fn deref(&self) -> &[T] {
        match self {
            Inline::Inline { len, items } => &items[..*len as usize],
            Inline::Heap(v) => v,
        }
    }
}

impl<T: Copy + Default, const N: usize> DerefMut for Inline<T, N> {
    fn deref_mut(&mut self) -> &mut [T] {
        match self {
            Inline::Inline { len, items } => &mut items[..*len as usize],
            Inline::Heap(v) => v,
        }
    }
}

impl<T: Copy + Default, const N: usize> From<Vec<T>> for Inline<T, N> {
    fn from(v: Vec<T>) -> Self {
        if v.len() <= N {
            Inline::from_slice(&v)
        } else {
            Inline::Heap(v)
        }
    }
}

impl<T: Copy + Default, const N: usize> From<&[T]> for Inline<T, N> {
    fn from(slice: &[T]) -> Self {
        Inline::from_slice(slice)
    }
}

impl<T: Copy + Default, const N: usize> From<&Vec<T>> for Inline<T, N> {
    fn from(v: &Vec<T>) -> Self {
        Inline::from_slice(v)
    }
}

impl<T: Copy + Default, const N: usize, const K: usize> From<[T; K]> for Inline<T, N> {
    fn from(array: [T; K]) -> Self {
        Inline::from_slice(&array)
    }
}

impl<T: Copy + Default, const N: usize> From<Inline<T, N>> for Vec<T> {
    fn from(list: Inline<T, N>) -> Vec<T> {
        match list {
            Inline::Inline { len, items } => items[..len as usize].to_vec(),
            Inline::Heap(v) => v,
        }
    }
}

impl<T: Copy + Default, const N: usize> FromIterator<T> for Inline<T, N> {
    fn from_iter<I: IntoIterator<Item = T>>(iter: I) -> Self {
        let mut list = Inline::new();
        for item in iter {
            list.push(item);
        }
        list
    }
}

impl<T: Copy + Default, const N: usize> Extend<T> for Inline<T, N> {
    fn extend<I: IntoIterator<Item = T>>(&mut self, iter: I) {
        for item in iter {
            self.push(item);
        }
    }
}

impl<'a, T: Copy + Default, const N: usize> IntoIterator for &'a Inline<T, N> {
    type Item = &'a T;
    type IntoIter = std::slice::Iter<'a, T>;
    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}

/// The owned iterator of an [`Inline`]: it yields copies, by position.
pub struct IntoIter<T: Copy + Default, const N: usize> {
    list: Inline<T, N>,
    next: usize,
}

impl<T: Copy + Default, const N: usize> Iterator for IntoIter<T, N> {
    type Item = T;
    fn next(&mut self) -> Option<T> {
        let item = self.list.get(self.next).copied();
        self.next += 1;
        item
    }
    fn size_hint(&self) -> (usize, Option<usize>) {
        let left = self.list.len().saturating_sub(self.next);
        (left, Some(left))
    }
}

impl<T: Copy + Default, const N: usize> ExactSizeIterator for IntoIter<T, N> {}

impl<T: Copy + Default, const N: usize> IntoIterator for Inline<T, N> {
    type Item = T;
    type IntoIter = IntoIter<T, N>;
    fn into_iter(self) -> IntoIter<T, N> {
        IntoIter {
            list: self,
            next: 0,
        }
    }
}

impl<T: Copy + Default + PartialEq, const N: usize> PartialEq for Inline<T, N> {
    fn eq(&self, other: &Self) -> bool {
        **self == **other
    }
}

impl<T: Copy + Default + Eq, const N: usize> Eq for Inline<T, N> {}

impl<T: Copy + Default + PartialEq, const N: usize> PartialEq<Vec<T>> for Inline<T, N> {
    fn eq(&self, other: &Vec<T>) -> bool {
        **self == **other
    }
}

impl<T: Copy + Default + PartialEq, const N: usize> PartialEq<[T]> for Inline<T, N> {
    fn eq(&self, other: &[T]) -> bool {
        **self == *other
    }
}

impl<T: Copy + Default + PartialEq, const N: usize, const K: usize> PartialEq<[T; K]>
    for Inline<T, N>
{
    fn eq(&self, other: &[T; K]) -> bool {
        **self == other[..]
    }
}

impl<T: Copy + Default + Hash, const N: usize> Hash for Inline<T, N> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        (**self).hash(state)
    }
}

impl<T: Copy + Default + fmt::Debug, const N: usize> fmt::Debug for Inline<T, N> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Debug::fmt(&**self, f)
    }
}

impl Instruction {
    pub fn gate(
        name: &str,
        qubits: impl Into<Inline<usize, 3>>,
        params: impl Into<Inline<f64, 3>>,
    ) -> Instruction {
        Instruction {
            name: Name::from(name),
            qubits: qubits.into(),
            params: params.into(),
            clbits: vec![],
        }
    }

    pub fn is_gate(&self) -> bool {
        !NON_GATES.contains(&self.name.as_str()) && self.clbits.is_empty()
    }

    /// The same instruction with qubit `i` replaced by `qubits[i]`.
    pub fn relabel(&self, qubits: &[usize]) -> Instruction {
        Instruction {
            name: self.name.clone(),
            qubits: self.qubits.iter().map(|&q| qubits[q]).collect(),
            params: self.params.clone(),
            clbits: self.clbits.clone(),
        }
    }
}

/// The cost the passes minimise: two-qubit gates first, then all gates.
pub fn cost(instructions: &[Instruction]) -> (usize, usize) {
    (
        instructions.iter().filter(|i| i.qubits.len() > 1).count(),
        instructions.len(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn known_names_are_borrowed_and_others_owned() {
        assert!(matches!(Name::from("cx").0, Cow::Borrowed(_)));
        assert!(matches!(
            Name::from(String::from("rzz")).0,
            Cow::Borrowed(_)
        ));
        let custom = Name::from("my_gate");
        assert!(matches!(custom.0, Cow::Owned(_)));
        assert_eq!(custom, "my_gate");
        assert_eq!(custom.to_string(), "my_gate");
        assert_eq!(format!("{:?}", Name::from("h")), "\"h\"");
        for &name in KNOWN_NAMES {
            assert!(matches!(Name::from(name).0, Cow::Borrowed(b) if b == name));
        }
    }

    #[test]
    fn short_lists_stay_inline_and_long_ones_spill() {
        let mut list: Inline<usize, 3> = Inline::new();
        assert!(list.is_empty());
        for q in 0..3 {
            list.push(q);
            assert!(matches!(list, Inline::Inline { .. }));
        }
        list.push(3);
        assert!(matches!(list, Inline::Heap(_)));
        assert_eq!(list, vec![0, 1, 2, 3]);
        let wide: Inline<usize, 3> = (0..7).collect();
        assert_eq!(&*wide, &[0, 1, 2, 3, 4, 5, 6]);
        assert_eq!(wide.clone().into_iter().sum::<usize>(), 21);
        let pair: Inline<usize, 3> = [4, 2].into();
        assert_eq!(Vec::from(pair.clone()), vec![4, 2]);
        assert_eq!(format!("{pair:?}"), "[4, 2]");
        // A list hashes and compares as its slice, wherever it is stored.
        let spilled: Inline<usize, 3> = Inline::Heap(vec![4, 2]);
        assert_eq!(pair, spilled);
        let hash = |list: &Inline<usize, 3>| {
            let mut h = std::collections::hash_map::DefaultHasher::new();
            list.hash(&mut h);
            h.finish()
        };
        assert_eq!(hash(&pair), hash(&spilled));
    }
}
