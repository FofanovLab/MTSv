//! Run quality control and deduplication processes on a batch of FASTQ files, producing a FASTA
//! file.

use bio::io::fastq::Reader;
use cue::pipeline;

use error::MtsvResult;
use itertools::Itertools;
use prep_config::{PrepConfig, TrimType};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::process::exit;

/// Execute QC processes on a given configuration, in parallel as much as possible.
pub fn run_prep(config: &PrepConfig) -> MtsvResult<()> {
    let mut processed: HashMap<Vec<u8>, Vec<usize>> = HashMap::new();
    let mut writer = BufWriter::new(try!(File::create(&config.outfile)));

    // chain together all of the FASTQ iterators
    let mut readers = Vec::new();
    for (i, &(ref path, _)) in config.infiles.iter().enumerate() {
        readers.push((i, try!(Reader::from_file(path))));
    }

    let reads = readers.into_iter()
        .flat_map(|(i, reader)| reader.records().map(move |r| (i, r)));

    // run the pipeline
    pipeline("prep reads",
             config.num_threads,
             reads,
             |(i, r)| {
        let r = match r {
            Ok(r) => r,
            Err(why) => {
                error!("Unable to read FASTQ file ({}): {:?}",
                       config.infiles[i].0.display(),
                       why);
                exit(13);
            },
        };

        // get the subsequence(s) from the QC processes
        let subseqs = process_read(r.seq(), r.qual(), i, config);

        (i, subseqs)
    },
             |(i, s)| {
        // deduplicate subsequences, keeping count of which file they came from
        for r in s {
            let counts = processed.entry(r)
                .or_insert_with(|| {
                    let mut v = Vec::with_capacity(config.infiles.len());
                    v.resize(config.infiles.len(), 0);
                    v
                });
            counts[i] += 1;
        }
    });

    let mut processed = processed.into_iter().collect::<Vec<_>>();
    processed.sort();

    // all reads should now be in order, so write out their results
    for (num, (mut read, counts)) in processed.into_iter().enumerate() {
        try!(write!(&mut writer, ">R{}", num + 1));
        for c in counts {
            try!(write!(&mut writer, "_{}", c));
        }
        try!(write!(&mut writer, "\n"));

        try!(writer.write_all(&mut read));
        try!(write!(&mut writer, "\n"));
    }

    Ok(())
}

/// Run quality control processes on a read.
///
/// * Trim adapters if enabled
/// * Trim sequence to uniform length (discard if removing adapters left it too short), segmenting
/// into multiple subsequences if enabled
/// * Filter out any low-quality subsequences if enabled (depending on minimum quality and # of low
/// quality bases tolerated)
pub fn process_read(seq: &[u8],
                    qual: &[u8],
                    file_index: usize,
                    config: &PrepConfig)
                    -> Vec<Vec<u8>> {

    let new_start = match config.adapters {
        Some(ref adapters) => {
            let tolerance = config.adapter_tolerance.unwrap();
            first_non_adapter_char(seq, adapters, tolerance)
        },
        None => 0,
    };

    let seq = &seq[new_start..];
    let qual = &qual[new_start..];

    let trim_len = match config.trim {
        TrimType::LcdFirstN(l) => l,
        TrimType::LcdQuality(l) => l,
        TrimType::Segment(l) => l,
    };

    let mut subseqs = Vec::new();
    if trim_len <= seq.len() {
        match config.trim {
            TrimType::LcdFirstN(n) => subseqs.push(lcd_trim(seq, qual, n)),
            TrimType::LcdQuality(n) => subseqs.push(lcd_q_trim(seq, qual, n)),
            TrimType::Segment(n) => subseqs.extend(segment_trim(seq, qual, n)),
        }
    }

    subseqs.into_iter()
        .filter(|&(_, qual)| match (config.min_quality, config.quality_threshold) {
            (Some(q), Some(n)) => {
                let q = config.infiles[file_index].1.qual_encoding.offset(q);

                is_high_enough_quality(qual, q, n)
            },
            _ => true,
        })
        .map(|(s, _)| Vec::from(s))
        .collect::<Vec<_>>()
}

/// Generates a 2-tuple of the first n base-pairs of the given sequence and quality.
pub fn lcd_trim<'a, 'b>(sequence: &'a [u8],
                        quality: &'b [u8],
                        trim_length: usize)
                        -> (&'a [u8], &'b [u8]) {
    assert!(sequence.len() == quality.len());

    (&sequence[0..trim_length], &quality[0..trim_length])
}

/// Given a sequence/quality pair, and a desired length, split the read into equal length segments,
/// discarding any tail portion which is shorter than the desired segment length.
pub fn segment_trim<'a, 'b>(sequence: &'a [u8],
                            quality: &'b [u8],
                            trim_length: usize)
                            -> Vec<(&'a [u8], &'b [u8])> {
    assert!(sequence.len() == quality.len());

    (0..sequence.len())
        .step(trim_length)
        .filter(|l| l + trim_length <= sequence.len())
        .map(|l| (&sequence[l..l + trim_length], &quality[l..l + trim_length]))
        .collect()
}

/// Generates a 2-tuple of the highest quality n-length subsequence and its quality scores.
pub fn lcd_q_trim<'a, 'b>(sequence: &'a [u8],
                          quality: &'b [u8],
                          trim_length: usize)
                          -> (&'a [u8], &'b [u8]) {
    assert!(sequence.len() == quality.len());

    let start = highest_q_start(quality, trim_length);
    let end = start + trim_length;

    (&sequence[start..end], &quality[start..end])
}

/// Find the start of the highest sum n-length subsequence of quality scores.
pub fn highest_q_start(quality: &[u8], length: usize) -> usize {
    assert!(0 < length && length <= quality.len());

    let mut max_start = 0;
    let mut max_sum: usize = (&quality[0..length]).iter().fold(0, |acc, &q| acc + q as usize);

    let mut curr_sum = max_sum;

    println!("quals: {:?}", quality);
    println!("length {} quality sum at start {}: {}", length, 0, curr_sum);
    for start in 1..(quality.len() - (length - 1)) {
        curr_sum -= quality[start] as usize;
        curr_sum += quality[start + (length - 1)] as usize;

        println!("length {} quality sum at start {}: {}",
                 length,
                 start,
                 curr_sum);

        if max_sum < curr_sum {
            max_sum = curr_sum;
            max_start = start;

            println!("MAX len {} quality sum at start {}: {}",
                     length,
                     start,
                     max_sum);
        }
    }

    println!("Found final max sum of {} at start {}", max_sum, max_start);

    max_start
}

/// Takes all substrings at the left of the sequences which are greater/equal in length than the
/// tolerance, and checks them against the provided set of adapters (which should include adapter
/// reverse complements, and their subsequences).
pub fn first_non_adapter_char(sequence: &[u8],
                              adapters: &HashSet<Vec<u8>>,
                              tolerance: usize)
                              -> usize {
    assert!(sequence.len() > tolerance);

    if adapters.len() == 0 {
        return 0;
    }

    for end in (tolerance..(sequence.len() + 1)).rev() {
        if adapters.contains(&sequence[..end]) {
            return end;
        }
    }

    0
}

/// Determines if a given list of base pair quality scores meets the configured thresholds.
pub fn is_high_enough_quality(quality: &[u8], min_quality: u8, tolerance: usize) -> bool {
    if min_quality == 0 {
        return true;
    }

    quality.iter().filter(|&&q| q < min_quality).count() <= tolerance
}

#[cfg(test)]
mod tests {
    use error::MtsvResult;
    use itertools::Itertools;
    use mktemp::Temp;
    use prep_config::*;
    use std::cmp;
    use std::path::Path;
    use super::*;

    #[test]
    fn test_prep_integration() {
        let base_dir = Path::new(env!("CARGO_MANIFEST_DIR")).to_path_buf();

        let outfile = Temp::new_file().unwrap();
        let outfile = outfile.to_path_buf();

        let adapter_path = base_dir.join("tests/prep/adapter.list");
        let adapters = read_adapters(&adapter_path, 10)
            .expect(&format!("Can't read adapters ({}).", adapter_path.display()));

        let expected_fasta_path = base_dir.join("tests/prep/expected_reads.fasta");
        let expected_fasta = file_to_bytes(&expected_fasta_path)
            .expect(&format!("Can't read expected FASTA ({:?})",
                             expected_fasta_path.display()));

        let files = vec![
            (base_dir.join("tests/prep/sample1.fastq"),
                FastqMetadata {
                read_len: 100,
                num_reads: 61927,
                qual_encoding: QualityEncoding::Sanger,
            }),
            (base_dir.join("tests/prep/sample3.fastq"),
            FastqMetadata {
                read_len: 100,
                num_reads: 52513,
                qual_encoding: QualityEncoding::Sanger,
            }),
            (base_dir.join("tests/prep/sample2.fastq"),
            FastqMetadata {
                read_len: 100,
                num_reads: 50402,
                qual_encoding: QualityEncoding::Sanger,
            }),
        ];

        let mut infiles = Vec::new();

        for (f, md) in files {
            let md_read = read_fastq_metadata(&f).unwrap();

            assert_eq!(md, md_read);

            infiles.push((f, md_read));
        }

        let config = PrepConfig {
            trim: TrimType::Segment(50),
            min_quality: Some(25),
            quality_threshold: Some(6),
            adapters: Some(adapters),
            adapter_tolerance: Some(10),
            num_threads: 4,
            infiles: infiles,
            outfile: outfile,
        };

        run_prep(&config).unwrap();

        let prepped_bytes = file_to_bytes(&config.outfile)
            .expect(&format!("Can't read prepared file ({:?})", config.outfile.display()));

        assert_eq!(expected_fasta.len(), prepped_bytes.len());

        let step = 100;
        for i in (0..expected_fasta.len()).step(step) {
            let end = cmp::min(i + step, expected_fasta.len());

            let expected = String::from_utf8_lossy(&expected_fasta[i..end]);
            let found = String::from_utf8_lossy(&prepped_bytes[i..end]);

            assert_eq!((i, expected), (i, found));
        }
    }

    fn file_to_bytes(p: &Path) -> MtsvResult<Vec<u8>> {
        use std::fs::File;
        use std::io::Read;

        let mut f = try!(File::open(p));
        let mut buf = Vec::new();
        try!(f.read_to_end(&mut buf));

        Ok(buf)
    }

    #[test]
    fn highest_q_same_length() {
        let test_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
        let start = highest_q_start(&test_list, test_list.len());

        assert_eq!(start, 0);
    }

    #[test]
    #[should_panic]
    fn test_highest_q_longer_desired_than_avail() {
        let test_list = [1, 2, 3, 4];
        highest_q_start(&test_list, 5);
    }

    #[test]
    #[should_panic]
    fn test_highest_q_zero_goal() {
        let test_list = [1, 2, 3, 4, 5];
        highest_q_start(&test_list, 0);
    }

    #[test]
    fn highest_at_start() {
        let test_list = [10, 9, 8, 7, 6, 5, 4, 3];
        assert_eq!(highest_q_start(&test_list, 8), 0);
        assert_eq!(highest_q_start(&test_list, 7), 0);
        assert_eq!(highest_q_start(&test_list, 6), 0);
        assert_eq!(highest_q_start(&test_list, 5), 0);
        assert_eq!(highest_q_start(&test_list, 4), 0);
        assert_eq!(highest_q_start(&test_list, 3), 0);
        assert_eq!(highest_q_start(&test_list, 2), 0);
    }

    #[test]
    fn highest_at_end() {
        let test_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
        assert_eq!(highest_q_start(&test_list, 10), 0);
        assert_eq!(highest_q_start(&test_list, 9), 1);
        assert_eq!(highest_q_start(&test_list, 8), 2);
        assert_eq!(highest_q_start(&test_list, 7), 3);
        assert_eq!(highest_q_start(&test_list, 6), 4);
        assert_eq!(highest_q_start(&test_list, 5), 5);
        assert_eq!(highest_q_start(&test_list, 4), 6);
        assert_eq!(highest_q_start(&test_list, 3), 7);
    }

    #[test]
    fn highest_in_middle() {
        let test_list = [1, 2, 3, 4, 5, 5, 4, 3, 2, 1];
        assert_eq!(highest_q_start(&test_list, 9), 0);
        assert_eq!(highest_q_start(&test_list, 8), 0);
        assert_eq!(highest_q_start(&test_list, 7), 1);
        assert_eq!(highest_q_start(&test_list, 6), 1);
        assert_eq!(highest_q_start(&test_list, 5), 2);
        assert_eq!(highest_q_start(&test_list, 4), 2);
        assert_eq!(highest_q_start(&test_list, 3), 3);
    }

    const ADAPTER: &'static [u8] = b"GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG";
    const TOLERANCE: usize = 10;

    #[test]
    fn adapters_whole_sequence() {
        let adapters = subadapters(&vec![ADAPTER], TOLERANCE);

        assert_eq!(ADAPTER.len(),
                   first_non_adapter_char(ADAPTER, &adapters, TOLERANCE));
    }

    #[test]
    fn adapters_none() {
        let adapters = subadapters(&vec![ADAPTER], TOLERANCE);

        let test_sequence = b"GAAGAGCCGTGCTTGGAAGAGCCGTGCTTGGAAGAGCCGTGCTTGGAAGAGCCGTGCTTG";
        assert_eq!(first_non_adapter_char(test_sequence, &adapters, TOLERANCE),
                   0);
    }

    #[test]
    fn adapters_empty() {
        let adapters = subadapters(&[], 1);
        let test_sequence = b"GAAGAGCCGTGCTTGGAAGAGCCGTGCTTGGAAGAGCCGTGCTTGGAAGAGCCGTGCTTG";
        assert_eq!(first_non_adapter_char(test_sequence, &adapters, TOLERANCE),
                   0);
    }

    #[test]
    fn adapters_present_but_below_tolerance() {
        let adapters = subadapters(&vec![ADAPTER], TOLERANCE);

        let test_sequence = b"TTCTGCTTGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        assert_eq!(first_non_adapter_char(test_sequence, &adapters, TOLERANCE),
                   0);
    }

    #[test]
    fn adapters_present_and_above_tolerance() {
        let adapters = subadapters(&vec![ADAPTER], TOLERANCE);

        let test_sequence = b"CTTCTGCTTGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        assert_eq!(first_non_adapter_char(test_sequence, &adapters, TOLERANCE),
                   10);

        let test_sequence = b"TCTTCCGATCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        assert_eq!(first_non_adapter_char(test_sequence, &adapters, TOLERANCE),
                   10);
    }

    #[test]
    fn quality_filter() {
        let qualities = [25, 25, 25, 25, 25, 24, 24, 24];

        assert!(!is_high_enough_quality(&qualities, 25, 2));
        assert!(!is_high_enough_quality(&qualities, 26, 5));

        assert!(is_high_enough_quality(&qualities, 24, 1));
        assert!(is_high_enough_quality(&qualities, 26, 8));

        assert!(is_high_enough_quality(&qualities, 0, 0));
    }

    #[test]
    fn lcd_trim_same_length() {
        let test_sequence = [b'A'; 10];
        let test_quality = [30; 10];

        let (seq, qual) = lcd_trim(&test_sequence, &test_quality, 10);

        assert_eq!(&test_sequence, seq);
        assert_eq!(&test_quality, qual);
    }

    #[test]
    fn lcd_trim_diff_length() {
        let test_seq = [b'A'; 10];
        let test_qual = [30; 10];

        let (seq, qual) = lcd_trim(&test_seq, &test_qual, 9);

        assert_eq!(seq, &[b'A'; 9]);
        assert_eq!(qual, &[30; 9]);

        let test_sequence = b"ABABABABABABABABABAB";
        let test_quality = [30; 20];

        let (seq, qual) = lcd_trim(test_sequence, &test_quality, 19);

        assert_eq!(seq, b"ABABABABABABABABABA");
        assert_eq!(qual, &[30; 19]);
    }

    #[test]
    fn lcdq_trim_same_length() {
        let test_sequence = b"AAAAABBBBBCCCCCDDDDD";
        let test_quality = [30, 30, 30, 30, 30, 31, 31, 31, 31, 31, 32, 32, 32, 32, 32, 33, 33,
                            33, 33, 33];

        let (result_seq, result_qual) = lcd_q_trim(test_sequence, &test_quality, 20);

        assert_eq!(result_seq, test_sequence);
        assert_eq!(result_qual, &test_quality);
    }

    #[test]
    fn lcdq_trim_diff_length() {
        let test_sequence = b"AAAAABBBBBCCCCCDDDDD";
        let test_quality = [30, 30, 30, 30, 30, 31, 31, 31, 31, 31, 32, 32, 32, 32, 32, 33, 33,
                            33, 33, 33];

        let (result_seq, result_qual) = lcd_q_trim(test_sequence, &test_quality, 10);

        assert_eq!(result_seq, b"CCCCCDDDDD");
        assert_eq!(result_qual, [32, 32, 32, 32, 32, 33, 33, 33, 33, 33]);
    }

    #[test]
    fn segment_trim_same_len() {
        let test_sequence = b"AAAAABBBBBCCCCCDDDDD";
        let test_quality = [30, 30, 30, 30, 30, 31, 31, 31, 31, 31, 32, 32, 32, 32, 32, 33, 33,
                            33, 33, 33];

        let results = segment_trim(test_sequence, &test_quality, 20);

        assert_eq!(results.len(), 1);

        let (result_sequence, result_quality) = results[0];

        assert_eq!(test_sequence, result_sequence);
        assert_eq!(&test_quality, result_quality);
    }

    #[test]
    fn segment_trim_even_division() {
        let test_sequence = b"AAAAABBBBBCCCCCDDDDD";
        let test_quality = [30, 30, 30, 30, 30, 31, 31, 31, 31, 31, 32, 32, 32, 32, 32, 33, 33,
                            33, 33, 33];

        let segments = segment_trim(test_sequence, &test_quality, 5);

        assert_eq!(segments.len(), 4);

        let (result_sequence, result_quality) = segments[0];
        assert_eq!(result_sequence, b"AAAAA");
        assert_eq!(result_quality, &[30, 30, 30, 30, 30]);

        let (result_sequence, result_quality) = segments[1];
        assert_eq!(result_sequence, b"BBBBB");
        assert_eq!(result_quality, &[31, 31, 31, 31, 31]);

        let (result_sequence, result_quality) = segments[2];
        assert_eq!(result_sequence, b"CCCCC");
        assert_eq!(result_quality, &[32, 32, 32, 32, 32]);

        let (result_sequence, result_quality) = segments[3];
        assert_eq!(result_sequence, b"DDDDD");
        assert_eq!(result_quality, &[33, 33, 33, 33, 33]);
    }

    #[test]
    fn segment_trim_uneven_division() {
        let test_sequence = b"AAAAABBBBBCCCCCDDDDD";
        let test_quality = [30, 30, 30, 30, 30, 31, 31, 31, 31, 31, 32, 32, 32, 32, 32, 33, 33,
                            33, 33, 33];

        let segments = segment_trim(test_sequence, &test_quality, 7);

        assert_eq!(segments.len(), 2);

        let (result_sequence, result_quality) = segments[0];
        assert_eq!(result_sequence, b"AAAAABB");
        assert_eq!(result_quality, [30, 30, 30, 30, 30, 31, 31]);

        let (result_sequence, result_quality) = segments[1];
        assert_eq!(result_sequence, b"BBBCCCC");
        assert_eq!(result_quality, [31, 31, 31, 32, 32, 32, 32]);

        let segments = segment_trim(test_sequence, &test_quality, 6);

        assert_eq!(segments.len(), 3);

        let (result_sequence, result_quality) = segments[0];
        assert_eq!(result_sequence, b"AAAAAB");
        assert_eq!(result_quality, [30, 30, 30, 30, 30, 31]);

        let (result_sequence, result_quality) = segments[1];
        assert_eq!(result_sequence, b"BBBBCC");
        assert_eq!(result_quality, [31, 31, 31, 31, 32, 32]);

        let (result_sequence, result_quality) = segments[2];
        assert_eq!(result_sequence, b"CCCDDD");
        assert_eq!(result_quality, [32, 32, 32, 33, 33, 33]);

        let segments = segment_trim(test_sequence, &test_quality, 3);
        assert_eq!(segments.len(), 6);

        let (result_sequence, result_quality) = segments[0];
        assert_eq!(result_sequence, b"AAA");
        assert_eq!(result_quality, [30, 30, 30]);

        let (result_sequence, result_quality) = segments[1];
        assert_eq!(result_sequence, b"AAB");
        assert_eq!(result_quality, [30, 30, 31]);

        let (result_sequence, result_quality) = segments[2];
        assert_eq!(result_sequence, b"BBB");
        assert_eq!(result_quality, [31, 31, 31]);

        let (result_sequence, result_quality) = segments[3];
        assert_eq!(result_sequence, b"BCC");
        assert_eq!(result_quality, [31, 32, 32]);

        let (result_sequence, result_quality) = segments[4];
        assert_eq!(result_sequence, b"CCC");
        assert_eq!(result_quality, [32, 32, 32]);

        let (result_sequence, result_quality) = segments[5];
        assert_eq!(result_sequence, b"DDD");
        assert_eq!(result_quality, [33, 33, 33]);
    }
}
