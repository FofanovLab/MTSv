//! The core metagenomic index used for queries.

use align::Aligner;

use bio::alphabets;
use bio::data_structures::bwt::bwt;
use bio::data_structures::fmindex::FMIndex;
use bio::data_structures::suffix_array::{SuffixArray, suffix_array};
use itertools::Itertools;
use ssw::{IDENT_W_PENALTY_NO_N_MATCH, Profile};
use std::cmp;
use std::collections::BTreeMap;
use std::fmt;
use std::fmt::{Debug, Formatter};
use std::hash::{Hash, Hasher};
use std::collections::hash_map::DefaultHasher;
use std::num::ParseIntError;
use std::str;
use std::u32;

/// Tuple struct to ensure GI/accession numbers don't get accidentally handled as tax IDs.
#[derive(Clone, Copy, Eq, PartialEq, PartialOrd, Ord, Hash, RustcEncodable, RustcDecodable, Debug)]
pub struct TaxId(pub u32);

/// Tuple struct to ensure taxonomic IDs don't get accidentally handled as GI/accession numbers.
#[derive(Clone, Copy, Eq, PartialEq, PartialOrd, Ord, Hash, RustcEncodable, RustcDecodable, Debug)]
pub struct Gi(pub u32);

/// Metadata about a region of the index, corresponding to a single sequence/GI/accession in the
/// original FASTA database file.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq, PartialOrd, Ord, RustcDecodable, RustcEncodable)]
struct Bin {
    gi: Gi,
    tax_id: TaxId,
    start: usize,
    end: usize,
}

/// Metagenomic index comprised of reference sequences concatenated together, an FM Index over the
/// concatenated sequences, and the metadata Bins to allow mapping absolute sequence offsets back
/// to GI/accession numbers and taxonomic IDs.
#[derive(Eq, Hash, PartialEq, RustcEncodable, RustcDecodable)]
pub struct MGIndex {
    fmindex: FMIndex,
    sequences: Sequence,
    bins: Vec<Bin>,
    suffix_array: SuffixArray,
}

impl Debug for MGIndex {
    fn fmt(&self, f: &mut Formatter) -> Result<(), fmt::Error> {
        let mut hasher = DefaultHasher::new();

        self.hash(&mut hasher);

        write!(f, "MGIndex {{ id: {}}}", hasher.finish())
    }
}

impl str::FromStr for TaxId {
    type Err = ParseIntError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match u32::from_str(s) {
            Ok(n) => Ok(TaxId(n)),
            Err(why) => Err(why),
        }
    }
}

impl str::FromStr for Gi {
    type Err = ParseIntError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match u32::from_str(s) {
            Ok(n) => Ok(Gi(n)),
            Err(why) => Err(why),
        }
    }
}

#[allow(missing_docs)]
pub type Sequence = Vec<u8>;

#[allow(missing_docs)]
pub type Database = BTreeMap<TaxId, Vec<(Gi, Sequence)>>;

/// The location within the index where a seed exact match was found.
#[derive(Copy, Clone, Debug, Eq, PartialEq, PartialOrd, Ord)]
struct SeedHit {
    reference_offset: usize,
    query_offset: usize,
}

impl SeedHit {
    /// Find the candidate alignment region for this seed hit, based on the query offset, the
    /// length of the original read, the edit distance tolerance, and the current GI bounds.
    pub fn candidate_indices(&self,
                             bin: &Bin,
                             read_len: usize,
                             edit_distance: usize)
                             -> Option<(usize, usize)> {
        let site = self.reference_offset;
        let seed_offset = self.query_offset;

        // the start of any alignment candidate needs to allow for some insertions at the beginning
        // but can't be earlier than the start of the GI in which this seed hit
        let start_offset = seed_offset + edit_distance;
        let cand_start = if site - start_offset < bin.start || start_offset > site {
            bin.start
        } else {
            site - start_offset
        };

        // same as the cand_start comment, but for the end of the current GI
        let cand_end = site + (read_len - seed_offset) + edit_distance;
        let cand_end = if cand_end > bin.end {
            bin.end
        } else {
            cand_end
        };

        // if:
        // these got swapped somehow
        // or we'd align against something outside the bin
        // or the candidate would be too short anyway
        if cand_start > cand_end || cand_start < bin.start || cand_end > bin.end ||
           cand_end - cand_start < read_len - edit_distance {
            None
        } else {
            Some((cand_start, cand_end))
        }
    }
}

/// A region of the reference sequences against which we may perform approximate alignment. Gets
/// expanded by adding successive `SeedHit`s.
#[derive(Copy, Clone, Eq, PartialEq)]
struct ReferenceCandidate<'rf> {
    reference_start: usize,
    reference_end_excl: usize,
    bin: Bin,
    num_seeds: usize,
    index: &'rf MGIndex,
}

impl<'rf> Debug for ReferenceCandidate<'rf> {
    fn fmt(&self, f: &mut Formatter) -> Result<(), fmt::Error> {
        write!(f,
               "ReferenceCandidate {{ start: {}, end: {}, bin: {:?}, num_seeds: {}, index: {:?}}}",
               self.reference_start,
               self.reference_end_excl,
               self.bin,
               self.num_seeds,
               self.index)
    }
}

impl<'rf> ReferenceCandidate<'rf> {
    /// Initialize a reference candidate with its first seed hit.
    fn new(seed_hit: SeedHit,
           bin: Bin,
           index: &'rf MGIndex,
           read_len: usize,
           edit_distance: usize)
           -> Option<Self> {

        let (ref_start, ref_end_excl) =
            match seed_hit.candidate_indices(&bin, read_len, edit_distance) {
                Some(r) => r,
                None => return None,
            };

        Some(ReferenceCandidate {
            reference_start: ref_start,
            reference_end_excl: ref_end_excl,
            bin: bin,
            num_seeds: 1,
            index: index,
        })
    }

    /// Returns a reference to the underlying candidate reference sequence.
    fn candidate_seq(&self) -> &'rf [u8] {
        &self.index.sequences[self.reference_start..self.reference_end_excl]
    }

    /// Attempts to merge another seed hit into this reference region. Succeeds if a candidate
    /// region derived from the new seed overlaps with the existing reference region. Fails if it
    /// would non-candidate portions of the reference into this candidate.
    fn add_seed_hit(&mut self,
                    seed_hit: SeedHit,
                    bin: &Bin,
                    read_len: usize,
                    edit_distance: usize)
                    -> Result<(), ()> {


        let (ref_start, ref_end_excl) =
            match seed_hit.candidate_indices(&bin, read_len, edit_distance) {
                Some(r) => r,
                None => return Err(()),
            };

        // check to see if this is even in the same GI
        if *bin != self.bin {
            Err(())
        }
        // check to see if the candidates from the new seed hit overlaps with the current candidate
        else if (self.reference_start <= ref_start && ref_start < self.reference_end_excl) ||
                  (self.reference_start < ref_end_excl && ref_end_excl <= self.reference_end_excl) {

            // since there's overlap, combine the two
            self.reference_start = cmp::min(self.reference_start, ref_start);
            self.reference_end_excl = cmp::max(self.reference_end_excl, ref_end_excl);

            // we filter and rank reference candidates by the number of seeds
            self.num_seeds += 1;

            Ok(())
        }
        // they're in the same bin, but the candidates don't overlap, we won't merge
        else {
            Err(())
        }
    }
}

impl MGIndex {
    // TODO test this function
    /// Identify all taxonomic IDs in this index which match against the query sequence within the
    /// specified edit distance.
    ///
    /// Process:
    ///
    /// 1. Generate a series of "seeds" (gapped subsequences) from the query sequence. The size and
    /// spacing of these are determined by the `seed_length` and `seed_gap` arguments.
    /// 2. Using the FM Index, find all locations where these seeds are present.
    /// 3. Coalesce all of the `SeedHit`s into combined `ReferenceCandidate`s representing all of
    /// the regions of the reference sequences against which we'll align the query sequence.
    /// 4. Sort all of the `ReferenceCandidate`s by the number of seeds present (we want to align
    /// the most likely regions first, as that will enable us to skip more regions later).
    /// 5. Use a SIMD-accelerated Smith-Waterman algorithm to align each reference candidate whose
    /// corresponding taxonomic ID hasn't already been found. When the score is within a threshold,
    /// perform a final edit-distance alignment, recording the taxonomic ID as "found" if it's
    /// equal to or lesser than the `edit_distance` argument.
    /// 6. Return the list of matching taxonomic IDs.

    pub fn matching_tax_ids(&self,
                            sequence: &[u8],
                            edit_distance: usize,
                            seed_length: usize,
                            seed_gap: usize,
                            min_seeds: usize)
                            -> Vec<TaxId> {

        // we need to later compare for edit distance where N's won't match against reference N's
        let seq_no_n = sequence.iter()
            .map(|b| {
                match *b {
                    b'N' => b'.',
                    _ => *b,
                }
            })
            .collect::<Vec<u8>>();

        let seeds = (0..(sequence.len() + 1 - seed_length)) // get all seed start indices
            .step(seed_gap)                                 // skip over any in between seed gap
            .map(|i| (i, &sequence[i..i + seed_length]));   // create a reference into the query

        // there are a few seeds which are SO prevalent they'll blow up memory usage if we don't
        // filter them out. in practice they have little impact on quality of results
        const SEED_FREQUENCY_TO_SKIP: usize = 100_000;

        // find all of the reference regions which we'll align against
        let reference_candidates = {
            let mut bin_locations = Vec::new();
            for (offset, seed) in seeds {
                // find everywhere this seed occurs in the reference database
                let interval = self.fmindex.backward_search(seed.iter());

                // if this seed is SUPER common, just skip it
                if interval.upper - interval.lower > SEED_FREQUENCY_TO_SKIP {
                    continue;
                }

                // track a new SeedHit for each value in ther suffix array interval
                bin_locations.extend((interval.lower..interval.upper).map(|i| {
                    SeedHit {
                        reference_offset: self.suffix_array[i],
                        query_offset: offset,
                    }
                }));
            }

            // merge all of the seed hits into candidate regions we can align against
            let mut refs =
                self.coalesce_seed_sites(&mut bin_locations,
                                         min_seeds,
                                         sequence.len(),
                                         edit_distance);

            // sort in reverse by number of seeds -- check the most promising locations first
            refs.sort_by(|a, b| b.num_seeds.cmp(&a.num_seeds));
            refs
        };

        let mut matches = Vec::new();

        let mut aligner = Aligner::new();

        let profile = Profile::new(sequence, &IDENT_W_PENALTY_NO_N_MATCH);

        for candidate in reference_candidates {
            // see if we've already found this tax ID
            if let Some(_) = matches.iter().find(|&&t| t == candidate.bin.tax_id) {
                continue;
            }

            // see if there's a match in the search candidate
            // if there is, record the hit tax id and then advance to the next candidate

            let cand_seq = candidate.candidate_seq();

            let score = profile.align_score(cand_seq, 1, 1);

            // -1 for substitution, -1 for gap open, -1 for gap extend
            // means that we need to allow for a hit to the alignment score of up to 1.5x editdist
            if score as usize >= sequence.len() - (edit_distance * 2) {

                // the SW check is faster (w/ SIMD) than the min_edit_distance check, so if we're
                // within an acceptable tolerance, now do the expensive check
                let edits = aligner.min_edit_distance(&seq_no_n, cand_seq);

                if edits as usize <= edit_distance {
                    matches.push(candidate.bin.tax_id);
                }
            }
        }

        matches
    }

    // TODO test this function
    /// Combine a series of `SeedHit`s into a series of `ReferenceCandidate`s.
    fn coalesce_seed_sites(&self,
                           seed_hits: &mut [SeedHit],
                           min_seeds: usize,
                           read_len: usize,
                           edit_distance: usize)
                           -> Vec<ReferenceCandidate> {
        seed_hits.sort();

        let mut curr_cand: Option<ReferenceCandidate> = None;
        let mut candidates = Vec::new();

        let mut bin_iter = self.bins.iter().peekable();
        // if there are no bins we have bigger problems
        let mut curr_bin = bin_iter.next().unwrap();

        for &mut sh in seed_hits {
            // if the site is ahead of the current bin, we need to advance the bin
            while curr_bin.end < sh.reference_offset {
                curr_bin = bin_iter.next().unwrap();
            }

            if let Some(mut cand) = curr_cand {

                if let Ok(()) = cand.add_seed_hit(sh, curr_bin, read_len, edit_distance) {
                    curr_cand = Some(cand);
                } else {
                    // if it wasn't added, it means that this seed hit is now past our current bin
                    // check if candidate has enough seeds, if so add to ref, set cand to None
                    if cand.num_seeds >= min_seeds {
                        candidates.push(cand);
                    }
                    curr_cand = None;
                }
            } else {
                curr_cand = ReferenceCandidate::new(sh, *curr_bin, self, read_len, edit_distance);
            }
        }

        candidates
    }

    /// Construct a new MGIndex from a series of reference sequences, concatenating all reference
    /// sequences and recording sequence boundaries and other metadata.
    pub fn new(reference: Database, sample_interval: u32) -> Self {
        info!("Concatenating all reference sequences and recording boundaries...");

        // concatenate all of the sequences, recording a new bin for each sequence
        let mut seq = Vec::new();
        let mut bins = Vec::new();
        for (tax_id, references) in reference {

            for (gi, reference) in references {
                let bin = Bin {
                    gi: gi,
                    tax_id: tax_id,
                    start: seq.len(),
                    end: seq.len() + reference.len(),
                };

                seq.extend_from_slice(&reference);
                bins.push(bin);
            }
        }

        // convert whole reference sequence to DNA5 alphabet
        for b in &mut seq {
            match *b {
                // skip capital N alphabet characters
                b'A' | b'C' | b'G' | b'T' | b'N' => (),
                b'a' => *b = b'A',
                b'c' => *b = b'C',
                b'g' => *b = b'G',
                b't' => *b = b'T',
                _ => *b = b'N',
            }
        }

        // suffix array requires a lexicographically smallest sentinel
        seq.push(b'$');
        seq.shrink_to_fit();
        info!("All reference sequences concatenated and boundaries recorded.");

        let alphabet = alphabets::dna::n_alphabet();

        info!("Building suffix array...");
        let suffix_array = suffix_array(&seq);
        info!("Suffix array constructed.");

        info!("Constructing Burrows-Wheeler Transform...");
        let bwt = bwt(&seq, &suffix_array);
        info!("BWT constructed.");

        info!("Constructing FM Index from SA and BWT...");
        let fmindex = FMIndex::new(bwt, sample_interval, &alphabet);
        info!("FMIndex constructed.");

        MGIndex {
            sequences: seq,
            fmindex: fmindex,
            bins: bins,
            suffix_array: suffix_array,
        }
    }
}

// this needs to be outside the test module so that integration tests can use it
#[cfg(test)]
pub fn random_database(num_taxa: u16,
                       num_gis: u16,
                       min_seq_size: usize,
                       max_seq_size: usize)
                       -> Database {
    use rand::{XorShiftRng, Rng};
    let mut rng = XorShiftRng::new_unseeded();

    let mut to_ret = BTreeMap::new();

    for _ in 0..num_taxa {
        let taxid = TaxId(rng.gen());
        let mut seqs = Vec::new();

        for _ in 0..num_gis {
            let gi = Gi(rng.gen());

            let mut seq = String::with_capacity(rng.gen_range(min_seq_size, max_seq_size));

            for _ in 0..seq.capacity() {
                let base = match rng.gen::<u8>() % 5 {
                    0 => 'A',
                    1 => 'C',
                    2 => 'G',
                    3 => 'T',
                    4 => 'N',
                    _ => unreachable!(),
                };
                seq.push(base);
            }

            seqs.push((gi, seq.into_bytes()));
        }

        to_ret.insert(taxid, seqs);
    }

    to_ret
}

#[cfg(test)]
mod test {
    use std::collections::BTreeMap;
    use super::*;
    use super::{Bin, ReferenceCandidate, SeedHit};

    #[test]
    #[should_panic]
    fn reference_candidate_non_overlapping() {
        let seed_hit = SeedHit {
            reference_offset: 110,
            query_offset: 1,
        };

        let seed_hit2 = SeedHit {
            reference_offset: 350,
            query_offset: 1,
        };

        let read_len = 50;
        let edits = 3;

        let db = random_database(10, 10, 500, 501);
        let index = MGIndex::new(db, 16);

        let bin = index.bins
            .iter()
            .filter(|b| b.start <= seed_hit.reference_offset && b.end > seed_hit.reference_offset)
            .next()
            .unwrap();

        let mut cand = ReferenceCandidate::new(seed_hit, *bin, &index, read_len, edits).unwrap();

        cand.add_seed_hit(seed_hit2, bin, read_len, edits).unwrap();
    }

    #[test]
    #[should_panic]
    fn reference_candidate_different_bin() {
        let seed_hit = SeedHit {
            reference_offset: 152,
            query_offset: 1,
        };

        let seed_hit2 = SeedHit {
            reference_offset: 350,
            query_offset: 1,
        };

        let read_len = 50;
        let edits = 3;

        let db = random_database(10, 10, 150, 151);
        let index = MGIndex::new(db, 16);

        if let Some(bin) = index.bins
            .iter()
            .filter(|b| b.start <= seed_hit.reference_offset && b.end > seed_hit.reference_offset)
            .next() {
            if let Some(bin2) = index.bins
                .iter()
                .filter(|b| {
                    b.start <= seed_hit2.reference_offset && b.end > seed_hit2.reference_offset
                })
                .next() {
                if let Some(mut cand) = ReferenceCandidate::new(seed_hit,
                                                                *bin,
                                                                &index,
                                                                read_len,
                                                                edits) {
                    // THIS is what should actually fail
                    cand.add_seed_hit(seed_hit2, bin2, read_len, edits).unwrap();
                }
            }
        }
    }

    #[test]
    fn reference_candidate() {
        let seed_hit = SeedHit {
            reference_offset: 110,
            query_offset: 1,
        };

        let read_len = 50;
        let edits = 3;

        let db = random_database(100, 200, 500, 1_000);
        let index = MGIndex::new(db, 16);

        let bin = index.bins
            .iter()
            .filter(|b| b.start <= seed_hit.reference_offset && b.end > seed_hit.reference_offset)
            .next()
            .unwrap();

        let mut cand = ReferenceCandidate::new(seed_hit, *bin, &index, read_len, edits).unwrap();

        let (expect_start, expect_end) = seed_hit.candidate_indices(bin, read_len, edits).unwrap();

        let found_seq = cand.candidate_seq();

        let found_ref_cand = ReferenceCandidate {
            reference_start: expect_start,
            reference_end_excl: expect_end,
            bin: *bin,
            num_seeds: 1,
            index: &index,
        };

        assert_eq!(found_ref_cand, cand);
        assert_eq!(found_seq, &index.sequences[expect_start..expect_end]);

        let seed_hit2 = SeedHit {
            reference_offset: 115,
            query_offset: 3,
        };

        cand.add_seed_hit(seed_hit2, bin, read_len, edits).unwrap();

        let (_, expect_end2) = seed_hit2.candidate_indices(bin, read_len, edits)
            .unwrap();

        assert_eq!(expect_start, cand.reference_start);
        assert_eq!(expect_end2, cand.reference_end_excl);
    }

    #[test]
    fn construct_index_lowercase() {
        let uppercase = random_database(100, 100, 150, 300);

        let lowercase: BTreeMap<_, _> = uppercase.iter()
            .map(|(taxon, seqs)| {
                let lc_seqs = seqs.iter()
                    .cloned()
                    .map(|(gi, seq)| {
                        (gi, String::from_utf8(seq).unwrap().to_lowercase().into_bytes())
                    })
                    .collect::<Vec<_>>();

                (*taxon, lc_seqs)
            })
            .collect();

        let uppercase = MGIndex::new(uppercase, 32);
        let lowercase = MGIndex::new(lowercase, 32);

        assert_eq!(uppercase, lowercase);
    }

    #[test]
    fn seed_hits_success() {
        let bin = Bin {
            gi: Gi(0),
            tax_id: TaxId(1),
            start: 100,
            end: 200,
        };

        let seed_hit = SeedHit {
            reference_offset: 110,
            query_offset: 1,
        };

        let read_len = 50;
        let edits = 3;
        let (cand_start, cand_end) = seed_hit.candidate_indices(&bin, read_len, edits).unwrap();

        assert!(cand_start < cand_end);
        assert!(cand_start >= bin.start);
        assert!(cand_end <= bin.end);
        assert!(cand_end - cand_start >= read_len + (2 * edits));

        let bin = Bin {
            gi: Gi(0),
            tax_id: TaxId(1),
            start: 100,
            end: 200,
        };

        let seed_hit = SeedHit {
            reference_offset: 180,
            query_offset: 25,
        };

        let read_len = 50;
        let edits = 3;
        let (cand_start, cand_end) = seed_hit.candidate_indices(&bin, read_len, edits).unwrap();

        assert!(cand_start < cand_end);
        assert!(cand_start >= bin.start);
        assert!(cand_end <= bin.end);
        assert!(cand_end - cand_start >= read_len - edits);
    }

    #[test]
    #[should_panic]
    fn seed_hits_fail() {
        let bin = Bin {
            gi: Gi(0),
            tax_id: TaxId(1),
            start: 100,
            end: 200,
        };

        let seed_hit = SeedHit {
            reference_offset: 90,
            query_offset: 1,
        };

        let read_len = 50;
        let edits = 3;
        let _ = seed_hit.candidate_indices(&bin, read_len, edits).unwrap();
    }
}
