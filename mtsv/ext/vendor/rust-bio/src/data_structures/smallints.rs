// Copyright 2014 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! A data structure for a sequence of small integers with a few big integers.
//! Small ints are stored in type S (e.g. a byte), big ints are stored separately (in type B) in a BTree.
//! The implementation provides vector-like operations on the data structure (e.g. retrieve a position,
//! add an integer, etc.).
//!
//! # Example
//!
//! ```
//! use bio::data_structures::smallints::SmallInts;
//! let mut smallints: SmallInts<u8, usize> = SmallInts::new();
//! smallints.push(3);
//! smallints.push(4);
//! smallints.push(255);
//! smallints.push(305093);
//! assert_eq!(smallints.get(0).unwrap(), 3);
//! smallints.set(0, 50000);
//! let values: Vec<usize> = smallints.iter().collect();
//! assert_eq!(values, [50000, 4, 255, 305093]);
//! ```

use std::iter::{Enumerate, repeat};
use std::mem::size_of;
use std::collections::BTreeMap;
use std::slice;

use num::{NumCast, Bounded, Integer, Num};
use num::traits::cast;


/// Data structure for storing a sequence of small integers with few big ones space efficiently
/// while supporting classical vector operations.
pub struct SmallInts<F: Integer + Bounded + NumCast + Copy, B: Integer + NumCast + Copy> {
    smallints: Vec<F>,
    bigints: BTreeMap<usize, B>,
}


impl<S: Integer + Bounded + NumCast + Copy, B: Integer + NumCast + Copy> SmallInts<S, B> {
    /// Create a new instance.
    pub fn new() -> Self {
        assert!(size_of::<S>() < size_of::<B>(),
                "S has to be smaller than B");
        SmallInts {
            smallints: Vec::new(),
            bigints: BTreeMap::new(),
        }
    }

    /// Create a new instance with a given capacity.
    pub fn with_capacity(n: usize) -> Self {
        assert!(size_of::<S>() < size_of::<B>(),
                "S has to be smaller than B");
        SmallInts {
            smallints: Vec::with_capacity(n),
            bigints: BTreeMap::new(),
        }
    }

    /// Create a new instance containing `n` times the integer `v` (and `v` is expected to be small).
    pub fn from_elem(v: S, n: usize) -> Self {
        assert!(size_of::<S>() < size_of::<B>(),
                "S has to be smaller than B");
        if v > cast(0).unwrap() {
            assert!(v < S::max_value(), "v has to be smaller than maximum value");
        }

        SmallInts {
            smallints: repeat(v).take(n).collect(),
            bigints: BTreeMap::new(),
        }
    }

    /// Return the integer at position `i`.
    pub fn get(&self, i: usize) -> Option<B> {
        if i < self.smallints.len() {
            self.real_value(i, self.smallints[i])
        } else {
            None
        }
    }

    /// Append `v` to the sequence. This will determine whether `v` is big or small and store it accordingly.
    pub fn push(&mut self, v: B) {
        let maxv: S = S::max_value();
        match cast(v) {
            Some(v) if v < maxv => self.smallints.push(v),
            _ => {
                let i = self.smallints.len();
                self.smallints.push(maxv);
                self.bigints.insert(i, v);
            }
        }
    }

    /// Set value of position `i` to `v`. This will determine whether `v` is big or small and store it accordingly.
    pub fn set(&mut self, i: usize, v: B) {
        let maxv: S = S::max_value();
        match cast(v) {
            Some(v) if v < maxv => self.smallints[i] = v,
            _ => {
                self.smallints[i] = maxv;
                self.bigints.insert(i, v);
            }
        }
    }

    /// Iterate over sequence. Values will be returned in the big integer type (`B`).
    pub fn iter(&self) -> Iter<S, B> {
        Iter {
            smallints: &self,
            items: self.smallints.iter().enumerate(),
        }
    }

    /// Decompress into a normal vector of big integers (type `B`).
    pub fn decompress(&self) -> Vec<B> {
        self.iter().collect()
    }

    /// Length of the sequence.
    pub fn len(&self) -> usize {
        self.smallints.len()
    }

    /// is the sequence empty?
    pub fn is_empty(&self) -> bool {
        self.smallints.is_empty()
    }

    fn real_value(&self, i: usize, v: S) -> Option<B> {
        if v < S::max_value() {
            cast(v)
        } else {
            self.bigints.get(&i).cloned()
        }
    }
}


/// Iterator over the elements of a SmallInts sequence.
pub struct Iter<'a, S, B>
    where S: 'a + Integer + Bounded + NumCast + Copy,
          B: 'a + Integer + NumCast + Copy,
          <S as Num>::FromStrRadixErr: 'a,
          <B as Num>::FromStrRadixErr: 'a
{
    smallints: &'a SmallInts<S, B>,
    items: Enumerate<slice::Iter<'a, S>>,
}


impl<'a, S, B> Iterator for Iter<'a, S, B>
    where S: 'a + Integer + Bounded + NumCast + Copy,
          B: 'a + Integer + NumCast + Copy,
          <S as Num>::FromStrRadixErr: 'a,
          <B as Num>::FromStrRadixErr: 'a
{
    type Item = B;

    fn next(&mut self) -> Option<B> {
        match self.items.next() {
            Some((i, &v)) => self.smallints.real_value(i, v),
            None => None,
        }
    }
}
