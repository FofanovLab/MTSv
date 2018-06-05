// Copyright 2014 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! FM-Index and FMD-Index for finding suffix array intervals matching a given pattern in linear time.

use std::iter::DoubleEndedIterator;

use data_structures::bwt::{Occ, Less, less, BWT};
use data_structures::suffix_array::SuffixArraySlice;
use alphabets::Alphabet;

/// The Fast Index in Minute space (FM-Index, Ferragina and Manzini, 2000) for finding suffix array
/// intervals matching a given pattern.
#[derive(Debug, Eq, Hash, PartialEq, RustcDecodable, RustcEncodable)]
pub struct FMIndex {
    bwt: BWT,
    less: Less,
    occ: Occ,
}

/// A suffix array interval.
#[derive(Debug, Copy, Clone)]
pub struct Interval {
    pub lower: usize,
    pub upper: usize,
}


impl Interval {
    /// Return the occurrence positions of the pattern as a slice of the suffix array.
    pub fn occ<'a>(&self, pos: &'a SuffixArraySlice) -> &'a [usize] {
        &pos[self.lower..self.upper]
    }
}


impl FMIndex {
    /// Construct a new instance of the FM index.
    ///
    /// # Arguments
    ///
    /// * `bwt` - the BWT
    /// * `k` - the sampling rate of the occ array: every k-th entry will be stored (higher k means
    ///   less memory usage, but worse performance)
    /// * `alphabet` - the alphabet of the underlying text, omitting the sentinel
    pub fn new(bwt: BWT, k: u32, alphabet: &Alphabet) -> Self {
        let less = less(&bwt, alphabet);
        let occ = Occ::new(&bwt, k);

        FMIndex {
            bwt: bwt,
            less: less,
            occ: occ,
        }
    }

    /// Perform backward search, yielding suffix array
    /// interval denoting exact occurences of the given pattern of length m in the text.
    /// Complexity: O(m).
    ///
    /// # Arguments
    ///
    /// * `pattern` - the pattern to search
    ///
    /// # Example
    ///
    /// ```
    /// use bio::data_structures::bwt::bwt;
    /// use bio::data_structures::fmindex::FMIndex;
    /// use bio::data_structures::suffix_array::suffix_array;
    /// use bio::alphabets::dna;
    ///
    /// let text = b"GCCTTAACATTATTACGCCTA$";
    /// let alphabet = dna::alphabet();
    /// let pos = suffix_array(text);
    /// let fm = FMIndex::new(bwt(text, &pos), 32, &alphabet);
    ///
    /// let pattern = b"TTA";
    /// let sai = fm.backward_search(pattern.iter());
    ///
    /// let occ = sai.occ(&pos);
    ///
    /// assert_eq!(occ, [3, 12, 9]);
    /// ```
    pub fn backward_search<'b, P: Iterator<Item = &'b u8> + DoubleEndedIterator>(&self,
                                                                                 pattern: P)
                                                                                 -> Interval {
        let (mut l, mut r) = (0, self.bwt.len() - 1);
        for &a in pattern.rev() {
            let less = self.less(a);
            l = less +
                if l > 0 {
                self.occ(l - 1, a)
            } else {
                0
            };
            r = less + self.occ(r, a) - 1;
        }

        Interval {
            lower: l,
            upper: r + 1,
        }
    }

    pub fn occ(&self, r: usize, a: u8) -> usize {
        self.occ.get(&self.bwt, r, a)
    }

    pub fn less(&self, a: u8) -> usize {
        self.less[a as usize]
    }

    /// Provide a reference to the underlying BWT.
    pub fn bwt(&self) -> &BWT {
        &self.bwt
    }
}
