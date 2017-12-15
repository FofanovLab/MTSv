//! A simple "minimum edit distance" sequence aligner with a reusable buffer.

use std::cmp::min;

/// An Aligner owns a buffer of data, and uses that to calculate the minimum edit distance with
/// which one sequence can be aligned against the other.
pub struct Aligner {
    buffer: Vec<u32>,
}

impl Aligner {
    /// Create a new Aligner. No allocations are performed until an alignment is requested.
    pub fn new() -> Self {
        Aligner { buffer: Vec::new() }
    }

    /// Find and return the minimum edit distance with which a needle can be aligned to a substring
    /// of haystack.
    ///
    /// Based on
    /// https://www.cs.jhu.edu/~langmea/resources/lecture_notes/variations_on_edit_dist.pdf.
    ///
    /// # Safety
    ///
    /// This method makes liberal use of `Vec::get_unchecked_mut`. All accesses are within bounds,
    /// but pay *very close* attention if modifying the indexing logic here.

    pub fn min_edit_distance(&mut self, p: &[u8], t: &[u8]) -> u32 {
        let dp_size = (p.len() + 1) * (t.len() + 1);
        let row_mult = t.len() + 1;

        let d = &mut self.buffer;

        d.resize(dp_size, 0);

        // Note: First row gets zeros.  First column initialized as usual in Levenshtein distance
        // calc.
        for i in 0..row_mult {
            unsafe {
                *d.get_unchecked_mut(i) = 0;
            }
        }

        // fill the first column with increasing numbers
        for row in 1..(p.len() + 1) {
            unsafe {
                *d.get_unchecked_mut(row * row_mult) = row as u32;
            }
        }

        for row in 1..(p.len() + 1) {
            for col in 1..(t.len() + 1) {

                unsafe {
                    // let needle_char = p[row - 1];
                    // let haystack_char = t[col - 1];
                    let needle_char = *p.get_unchecked(row - 1);
                    let haystack_char = *t.get_unchecked(col - 1);

                    // do the characters at this cell match? if not, potentially add 1 to edit dist
                    let delta = if needle_char != haystack_char { 1 } else { 0 };

                    // determine score weights for insertion, deletion, substitution
                    let diag = ((row - 1) * row_mult) + (col - 1);
                    let up = ((row - 1) * row_mult) + col;
                    let left = (row * row_mult) + (col - 1);

                    // d[(row * row_mult) + col] = min(d[diag] + delta, min(d[up] + 1, d[left] +
                    // 1));
                    let new_current = min(*d.get_unchecked(diag) + delta,
                                          min(*d.get_unchecked(up) + 1,
                                              *d.get_unchecked(left) + 1));

                    let current = d.get_unchecked_mut((row * row_mult) + col);

                    *current = new_current;
                }

            }
        }

        // get the minimum value in the last row
        let last_row = &d[(dp_size - (t.len() + 1))..dp_size];
        last_row.iter().map(|s| *s).min().unwrap()
    }
}

#[cfg(test)]
mod test {
    use super::Aligner;

    fn check_test(needle: &[u8], haystack: &[u8], expected_edits: u32) {
        let mut aligner = Aligner::new();

        let dist = aligner.min_edit_distance(needle, haystack);

        assert_eq!(dist, expected_edits);
    }

    #[test]
    fn test_from_jupyter_notebook() {
        let needle = b"TACGTCAGC";
        let haystack = b"AACCCTATGTCATGCCTTGGA";

        check_test(needle, haystack, 2);
    }

    #[test]
    fn test_exact_full() {
        let needle = b"ACGACTAGTTATAAAAATTCNACTCCANTTAGCTCCCTACTTTCCGAGAG";
        let haystack = b"ACGACTAGTTATAAAAATTCNACTCCANTTAGCTCCCTACTTTCCGAGAG";

        check_test(needle, haystack, 0);
    }

    #[test]
    fn test_exact_partial() {
        let needle = b"AAAAAT";
        let haystack = b"ACGACTAGTTATAAAAATTCNACTCCANTTAGCTCCCTACTTTCCGAGAG";

        check_test(needle, haystack, 0);
    }

    #[test]
    fn test_empty() {
        let needle = b"";
        let haystack = b"ACGACTAGTTATAAAAATTCNACTCCANTTAGCTCCCTACTTTCCGAGAG";

        check_test(needle, haystack, 0);
    }

    #[test]
    fn test_nomatches() {
        let needle = b"*********";
        let haystack = b"ACGACTAGTTATAAAAATTCNACTCCANTTAGCTCCCTACTTTCCGAGAG";

        check_test(needle, haystack, 9);
    }

    #[test]
    fn test_middle_edits() {
        let needle = b"ANNGTTCNGNT";
        let haystack = b"ACGACTAGTTATAAAAATTCNACTCCANTTAGCTCCCTACTTTCCGAGAG";

        check_test(needle, haystack, 5);
    }

    #[test]
    fn test_begin_edits() {
        let needle = b"***GTTATAA";
        let haystack = b"ACGACTAGTTATAAAAATTCNACTCCANTTAGCTCCCTACTTTCCGAGAG";

        check_test(needle, haystack, 3);
    }

    #[test]
    fn test_end_edits() {
        let needle = b"GTTATAA***";
        let haystack = b"ACGACTAGTTATAAAAATTCNACTCCANTTAGCTCCCTACTTTCCGAGAG";

        check_test(needle, haystack, 3);
    }
}
