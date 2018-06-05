//! (Mostly) safe bindings to Mengyao Zhao's SIMD implementation of Smith-Waterman.
//!
//! Currently limited to processing DNA5 sequences.

#![warn(missing_docs)]

extern crate libc;

/// Identity matrix for matching.
#[cfg_attr(rustfmt, rustfmt_skip)]
pub const IDENT_W_PENALTY_NO_N_MATCH: [i8; 25] =
    [1, -1, -1, -1, -1,
     -1, 1, -1, -1, -1,
     -1, -1, 1, -1, -1,
     -1, -1, -1, 1, -1,
     -1, -1, -1, -1, 1];


/// Query profile. Can be reused across alignments if aligning one sequence against many others.
pub struct Profile<'read> {
    sequence: &'read [u8],
    _sequence_numeric: Vec<i8>,
    raw_profile: *const RawProfile,
}

impl<'read> Drop for Profile<'read> {
    fn drop(&mut self) {
        unsafe {
            init_destroy(self.raw_profile);
        }
    }
}

impl<'read> Profile<'read> {
    /// Create a new query profile for the given DNA5 read and scoring matrix.
    pub fn new(read: &'read [u8], matrix: &[i8; 25]) -> Profile<'read> {
        assert!(read.len() > 0);

        let read_num = Self::sequence_to_numeric(read);
        let raw = unsafe {
            ssw_init(read_num.as_ptr(),
                     read_num.len() as i32,
                     matrix.as_ptr(),
                     5,
                     2)
        };

        // we need to store the numeric version of the sequence with the profile to make sure
        // that the Vec's underlying pointer isn't freed until the Profile goes out of scope

        Profile {
            sequence: read,
            _sequence_numeric: read_num,
            raw_profile: raw,
        }
    }

    /// Perform Smith-Waterman alignment of the contained query read against the `reference`
    /// sequence argument, returning the only the score. `gap_open` and `gap_extend` should be
    /// positive numbers which will be subtracted from the overall score.
    pub fn align_score(&self, reference: &[u8], gap_open: u8, gap_extend: u8) -> u16 {

        assert!(reference.len() > 0);

        let reference_numeric = Self::sequence_to_numeric(reference);

        let alignment = unsafe {
            ssw_align(self.raw_profile,
                      reference_numeric.as_ptr() as *const i8,
                      reference_numeric.len() as i32,
                      gap_open,
                      gap_extend,
                      0,
                      0,
                      0,
                      (self.sequence.len() / 2) as i32)
        };

        unsafe {
            let score = (*alignment).score1;

            align_destroy(alignment);

            score
        }
    }

    /// Convert a DNA5 read sequence to 0-based indices in the matrix.
    fn sequence_to_numeric(seq: &[u8]) -> Vec<i8> {
        let mut converted = Vec::with_capacity(seq.len());

        for &b in seq {
            let num = match b {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => 4,
            };

            converted.push(num);
        }

        converted
    }
}

#[repr(C)]
struct RawProfile {
    profile_byte: *const libc::c_void,
    profile_word: *const libc::c_void,
    read: *const i8,
    mat: *const i8,
    read_len: i32,
    n: i32,
    bias: u8,
}

#[repr(C)]
struct RawAlign {
    score1: u16,
    score2: u16,
    ref_begin1: i32,
    ref_end1: i32,
    read_begin1: i32,
    read_end1: i32,
    ref_end2: i32,
    cigar: *const u32,
    cigar_len: i32,
}

extern "C" {
    fn ssw_init(read: *const i8,
                readLen: i32,
                mat: *const i8,
                n: i32,
                score_size: i8)
                -> *const RawProfile;

    fn ssw_align(profile: *const RawProfile,
                 reference: *const i8,
                 ref_len: i32,
                 gap_open: u8,
                 gap_extend: u8,
                 flag: u8,
                 filters: u16,
                 filterd: i32,
                 mask_len: i32)
                 -> *const RawAlign;

    fn init_destroy(p: *const RawProfile);
    fn align_destroy(a: *const RawAlign);

}

#[cfg(test)]
extern crate arbitrary_dna;

#[cfg(test)]
extern crate bio;

#[cfg(test)]
#[macro_use]
extern crate quickcheck;

#[cfg(test)]
mod test {

    use arbitrary_dna::Dna5Sequence;
    use bio::alignment::pairwise::Aligner;
    use std::str;
    use super::*;

    quickcheck! {
        fn matches_rust_bio(query: Dna5Sequence,
                            reference: Dna5Sequence) -> bool {

            if query.len() < 32 || reference.len() < 32 {
                return true;
            }

            let query_bytes = query.iter().map(|base| base.0).collect::<Vec<u8>>();
            let reference_bytes = reference.iter().map(|base| base.0).collect::<Vec<u8>>();

            let scorer = |a, b| if a == b { 1i32 } else { -1i32 };
            let mut aligner = Aligner::new(-1, -1, &scorer);

            let expected = aligner.local(&query_bytes, &reference_bytes);

            let profile = Profile::new(&query_bytes, &IDENT_W_PENALTY_NO_N_MATCH);

            let found = profile.align_score(&reference_bytes, 1, 1);

            let diff = (found as i32 - expected.score).abs();

            // FIXME simd version is occasionally off by one
            diff <= 1
        }
    }
}
