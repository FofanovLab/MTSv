// Copyright 2014 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! The Burrows-Wheeler-Transform and related data structures.
//! The implementation is based on the lecture notes
//! "Algorithmen auf Sequenzen", Kopczynski, Marschall, Martin and Rahmann, 2008 - 2015.

use std::iter::repeat;

use utils::prescan;
use alphabets::Alphabet;
use data_structures::suffix_array::SuffixArraySlice;

pub type BWT = Vec<u8>;
pub type BWTSlice = [u8];
pub type Less = Vec<usize>;
pub type BWTFind = Vec<usize>;


/// Calculate Burrows-Wheeler-Transform of the given text of length n.
/// Complexity: O(n).
///
/// # Arguments
///
/// * `text` - the text ended by sentinel symbol (being lexicographically smallest)
/// * `pos` - the suffix array for the text
///
/// # Example
///
/// ```
/// use bio::data_structures::suffix_array::suffix_array;
/// use bio::data_structures::bwt::bwt;
/// let text = b"GCCTTAACATTATTACGCCTA$";
/// let pos = suffix_array(text);
/// let bwt = bwt(text, &pos);
/// assert_eq!(bwt, b"ATTATTCAGGACCC$CTTTCAA");
/// ```
pub fn bwt(text: &[u8], pos: &SuffixArraySlice) -> BWT {
    assert!(text.len() == pos.len());
    let n = text.len();
    let mut bwt: BWT = repeat(0).take(n).collect();
    for r in 0..n {
        let p = pos[r];
        bwt[r] = if p > 0 {
            text[p - 1]
        } else {
            text[n - 1]
        };
    }

    bwt
}


/// Calculate the inverse of a BWT of length n, which is the original text.
/// Complexity: O(n).
///
/// # Arguments
///
/// * `bwt` - the BWT
pub fn invert_bwt(bwt: &BWTSlice) -> Vec<u8> {
    let alphabet = Alphabet::new(bwt);
    let n = bwt.len();
    let bwtfind = bwtfind(bwt, &alphabet);
    let mut inverse = Vec::with_capacity(n);

    let mut r = bwtfind[0];
    for _ in 0..n {
        r = bwtfind[r];
        inverse.push(bwt[r]);
    }

    inverse
}

const OCC_BINS: usize = 6;

/// An occurence array implementation.
#[derive(Debug, Eq, Hash, PartialEq, RustcDecodable, RustcEncodable)]
pub struct Occ {
    occ: Vec<[usize; OCC_BINS]>,
    k: u32,
}

impl Occ {
    /// Calculate occ array with sampling from BWT of length n.
    /// Time complexity: O(n).
    /// Space complexity: O(n / k * A) with A being the alphabet size.
    /// Alphabet size is determined on the fly from the BWT.
    /// For large texts, it is therefore advisable to transform
    /// the text before calculating the BWT (see alphabets::rank_transform).
    ///
    /// # Arguments
    ///
    /// * `bwt` - the BWT
    /// * `k` - the sampling rate: every k-th entry will be stored
    pub fn new(bwt: &BWTSlice, k: u32) -> Self {
        let n = bwt.len();

        let mut occ = Vec::with_capacity(n / k as usize);
        let mut curr_occ = [0; OCC_BINS];

        for (i, &c) in bwt.iter().enumerate() {
            curr_occ[Occ::array_index(c)] += 1;
            if i % k as usize == 0 {
                occ.push(curr_occ.clone());
            }
        }

        Occ { occ: occ, k: k }
    }

    fn array_index(a: u8) -> usize {
        match a {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            b'N' => 4,
            _ => 5,
        }
    }

    /// Get occurrence count of symbol a in BWT[..r+1].
    /// Complexity: O(k).
    #[inline]
    pub fn get(&self, bwt: &BWTSlice, r: usize, a: u8) -> usize {
        let i = r / self.k as usize;

        // this is the count up to this point
        let checkpoint = self.occ[i][Occ::array_index(a)];

        // start counting in the BWT after the last sampled checkpoint
        let start = (i * self.k as usize) + 1;
        // end counting at the point of our query
        let end = r + 1;

        let mut count = 0u32;
        for &b in &bwt[start..end] {
            if b == a {
                count += 1;
            }
        }

        checkpoint + (count as usize)
    }
}


/// Calculate the less array for a given BWT. Complexity O(n).
pub fn less(bwt: &BWTSlice, alphabet: &Alphabet) -> Less {
    let m = alphabet.max_symbol() as usize + 2;
    let mut less: Less = repeat(0)
                             .take(m)
                             .collect();
    for &c in bwt.iter() {
        less[c as usize] += 1;
    }
    // calculate +-prescan
    prescan(&mut less[..], 0, |a, b| a + b);

    less
}


/// Calculate the bwtfind array needed for inverting the BWT. Complexity O(n).
pub fn bwtfind(bwt: &BWTSlice, alphabet: &Alphabet) -> BWTFind {
    let n = bwt.len();
    let mut less = less(bwt, alphabet);

    let mut bwtfind: BWTFind = repeat(0).take(n).collect();
    for (r, &c) in bwt.iter().enumerate() {
        bwtfind[less[c as usize]] = r;
        less[c as usize] += 1;
    }

    bwtfind
}


#[cfg(test)]
mod tests {
    use super::{bwtfind, bwt, invert_bwt, Occ};
    use data_structures::suffix_array::suffix_array;
    use alphabets::Alphabet;

    #[test]
    fn test_bwtfind() {
        let text = b"cabca$";
        let alphabet = Alphabet::new(b"abc$");
        let pos = suffix_array(text);
        let bwt = bwt(text, &pos);
        let bwtfind = bwtfind(&bwt, &alphabet);
        assert_eq!(bwtfind, vec![5, 0, 3, 4, 1, 2]);
    }

    #[test]
    fn test_invert_bwt() {
        let text = b"cabca$";
        let pos = suffix_array(text);
        let bwt = bwt(text, &pos);
        let inverse = invert_bwt(&bwt);
        assert_eq!(inverse, text);
    }

    #[test]
    fn test_occ() {
        let bwt = vec![b'N', b'C', b'C', b'N', b'A', b'G'];
        let occ = Occ::new(&bwt, 3);
        assert_eq!(occ.get(&bwt, 4, b'A'), 1);
        assert_eq!(occ.get(&bwt, 4, b'C'), 2);
    }
}
