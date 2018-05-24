// Copyright 2014 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Backward oracle matching algorithm.
//! Best-case complexity: O(n / m) with pattern of length m and text of length n.
//! Worst case complexity: O(n * m).
//!
//! # Example
//!
//! ```
//! use bio::pattern_matching::bom::BOM;
//! let text = b"ACGGCTAGGAAAAAGACTGAGGACTGAAAA";
//! let pattern = b"GAAAA";
//! let bom = BOM::new(pattern);
//! let occ: Vec<usize> = bom.find_all(text).collect();
//! assert_eq!(occ, [8, 25]);
//! ```


use std::iter::repeat;
use utils::TextSlice;

use vec_map::VecMap;


/// Backward oracle matching algorithm.
pub struct BOM {
    m: usize,
    table: Vec<VecMap<usize>>,
}


impl BOM {
    /// Create a new instance for a given pattern.
    pub fn new(pattern: TextSlice) -> Self {
        let m = pattern.len();
        let maxsym = *pattern.iter().max().expect("Expecting non-empty pattern.") as usize;
        let mut table: Vec<VecMap<usize>> = Vec::with_capacity(m);
        // init suffix table, initially all values unknown
        // suff[i] is the state in which the longest suffix of
        // pattern[..i+1] ends that does not end in i
        let mut suff: Vec<Option<usize>> = repeat(None).take(m + 1).collect();

        for (j, &b) in pattern.iter().rev().enumerate() {
            let i = j + 1;
            let a = b as usize;
            let mut delta = VecMap::with_capacity(maxsym);
            // reading symbol a leads into state i (this is an inner edge)
            delta.insert(a, i);
            // now, add edges for substrings ending with a
            let mut k = suff[i - 1];

            // for this iterate over the known suffixes until
            // reaching an edge labelled with a or the start
            while let Some(k_) = k {
                if table[k_].contains_key(a) {
                    break;
                }
                table[k_].insert(a, i);
                k = suff[k_];
            }

            // the longest suffix is either 0 or the state
            // reached by the edge labelled with a
            suff[i] = Some(match k {
                Some(k) => *table[k].get(a).unwrap(),
                None => 0,
            });

            table.push(delta);
        }

        BOM {
            m: m,
            table: table,
        }
    }

    fn delta(&self, q: usize, a: u8) -> Option<usize> {
        if q >= self.table.len() {
            None
        } else {
            match self.table[q].get(a as usize) {
                Some(&q) => Some(q),
                None => None,
            }
        }
    }

    /// Find all matches of the pattern in the given text. Matches are returned as an iterator over start positions.
    pub fn find_all<'a>(&'a self, text: TextSlice<'a>) -> Matches {
        Matches {
            bom: self,
            text: text,
            window: self.m,
        }
    }
}


/// Iterator over start positions of matches.
pub struct Matches<'a> {
    bom: &'a BOM,
    text: TextSlice<'a>,
    window: usize,
}


impl<'a> Iterator for Matches<'a> {
    type Item = usize;

    fn next(&mut self) -> Option<usize> {
        while self.window <= self.text.len() {
            let (mut q, mut j) = (Some(0), 1);
            while j <= self.bom.m {
                match q {
                    Some(q_) => {
                        q = self.bom.delta(q_, self.text[self.window - j]);
                        j += 1;
                    }
                    None => break,
                }
            }
            // putative start position
            let i = self.window - self.bom.m;
            self.window += self.bom.m + 2 - j;
            if q.is_some() {
                // return match
                return Some(i);
            }
        }
        None
    }
}


#[cfg(test)]
mod tests {
    use super::BOM;

    #[test]
    fn test_delta() {
        let pattern = b"qnnnannan"; // reverse of nannannnq
        let bom = BOM::new(pattern);
        assert_eq!(bom.delta(0, b'n'), Some(1));
        assert_eq!(bom.delta(1, b'a'), Some(2));
        assert_eq!(bom.delta(2, b'n'), Some(3));
        assert_eq!(bom.delta(3, b'n'), Some(4));
        assert_eq!(bom.delta(4, b'a'), Some(5));
        assert_eq!(bom.delta(5, b'n'), Some(6));
        assert_eq!(bom.delta(6, b'n'), Some(7));
        assert_eq!(bom.delta(7, b'n'), Some(8));
        assert_eq!(bom.delta(8, b'q'), Some(9));

        assert_eq!(bom.delta(0, b'a'), Some(2));
        assert_eq!(bom.delta(0, b'q'), Some(9));
        assert_eq!(bom.delta(1, b'n'), Some(4));
        assert_eq!(bom.delta(1, b'q'), Some(9));
        assert_eq!(bom.delta(4, b'n'), Some(8));
        assert_eq!(bom.delta(4, b'q'), Some(9));
        bom.delta(9, b'a');
    }
}
