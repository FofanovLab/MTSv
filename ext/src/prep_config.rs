//! Configuration for readprep.

use bio::alphabets::dna::RevComp;
use bio::io::fastq;
use bio::io::fastq::Reader;
use clap::{App, Arg, ArgGroup, ArgMatches};

use error::MtsvResult;
use std::cmp::min;
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::{Path, PathBuf};

/// The configuration for a particular readprep run.
#[derive(Debug, Eq, PartialEq, RustcEncodable)]
pub struct PrepConfig {
    /// How to trim sequences.
    pub trim: TrimType,
    /// The minimum quality to define, if filtering by quality.
    pub min_quality: Option<u8>,
    /// The number of low-quality bases to tolerate, if filtering by quality.
    pub quality_threshold: Option<usize>,
    /// The adapters to check for at the beginning of sequences, if filtering adapters.
    pub adapters: Option<HashSet<Vec<u8>>>,
    /// The number of adapter-matching bases to ignore at the beginning of a sequence, if filtering
    /// adapters.
    pub adapter_tolerance: Option<usize>,
    /// The number of threads to use for parallel execution.
    pub num_threads: usize,
    /// The input FASTQ files, along with the metadata read about them.
    pub infiles: Vec<(PathBuf, FastqMetadata)>,
    /// The path to write the FASTA results file to.
    pub outfile: PathBuf,
}

/// Which type of length-homogenization (trimming) to use on the reads in a particular file.
#[derive(Debug, Eq, PartialEq, RustcEncodable)]
pub enum TrimType {
    /// Just take the first N bases of the read sequence.
    LcdFirstN(usize),
    /// Find the highest quality, contiguous N bases of the read sequence, and use that.
    LcdQuality(usize),
    /// Split the read sequence into segments of N size, discarding any trailing bases.
    Segment(usize),
}

/// The inferred quality encoding used by a particular file.
#[derive(Debug, Eq, PartialEq, RustcEncodable)]
pub enum QualityEncoding {
    /// Illumina 1.5 (offset = 64)
    Illumina15,
    /// Illumina 1.3 (offset = 67)
    Illumina13,
    /// Solexa (offset 59)
    Solexa,
    /// Sanger (offset 33)
    Sanger,
}

impl QualityEncoding {
    /// Determine the encoding-specific value of an abstract quality value.
    pub fn offset(&self, quality: u8) -> u8 {
        match *self {
            QualityEncoding::Illumina13 => quality + 64,
            QualityEncoding::Illumina15 => quality + 67,
            QualityEncoding::Solexa => quality + 59,
            QualityEncoding::Sanger => quality + 33,
        }
    }
}

/// The metadata obtained about a particular FASTQ input file.
#[derive(Debug, Eq, PartialEq, RustcEncodable)]
pub struct FastqMetadata {
    /// The read length of all reads in the file. Assumes files have uniform read length.
    pub read_len: usize,
    /// The number of sequences present in the file.
    pub num_reads: usize,
    /// The inferred quality encoding.
    pub qual_encoding: QualityEncoding,
}

/// Run an initial parse of a file, determining it's read length, the number of reads contained,
/// and a best guess at what quality score encoding is used.
pub fn read_fastq_metadata(p: &Path) -> MtsvResult<FastqMetadata> {
    let rdr = try!(Reader::from_file(p));

    let mut count = 0;
    let mut read_len = 10_000_000;
    let mut min_quality = 255u8;

    for record in rdr.records() {
        let record = try!(record);

        count += 1;
        read_len = min(read_len, record.seq().len());
        let curr_min_qual = record.qual().iter().min().expect("invalid quality string found");

        min_quality = min(min_quality, *curr_min_qual);
    }

    let encoding = match min_quality {
        67...255 => QualityEncoding::Illumina15,
        64...66 => QualityEncoding::Illumina13,
        59...63 => QualityEncoding::Solexa,
        _ => QualityEncoding::Sanger,
    };

    Ok(FastqMetadata {
        read_len: read_len,
        num_reads: count,
        qual_encoding: encoding,
    })
}

/// Generate a list of adapters to check against from a given path.
pub fn read_adapters(f: &Path, tolerance: usize) -> MtsvResult<HashSet<Vec<u8>>> {
    let f = try!(File::open(f));
    let mut rdr = BufReader::new(f);
    let mut buf = String::new();

    try!(rdr.read_to_string(&mut buf));

    Ok(subadapters(&buf.lines().map(|l| l.trim().as_bytes()).collect::<Vec<_>>(),
                   tolerance))
}

/// Take a list of adapters, and generate a set of "subadapters" within a certain tolerance.
///
/// For example, if "ACGTACGTACGT" is an adapter, and we care about removing adapters that are at
/// least 4 bases in length, we would have the list:
///
/// * "ACGT"
/// * "TACGT"
/// * "GTACGT"
/// * "CGTACGT"
/// * "ACGTACGT"
/// * "TACGTACGT"
/// * "GTACGTACGT"
/// * "CGTACGTACGT"
/// * "ACGTACGTACGT"
///
/// Which would then be used for matching against the initial N bases of a sequence to determine
/// whether it had been contaminated by an adapter.
pub fn subadapters(full_adapters: &[&[u8]], tolerance: usize) -> HashSet<Vec<u8>> {
    let mut adapters = HashSet::new();
    let revcomp = RevComp::new();

    for a in full_adapters {
        let rev = revcomp.get(a);

        // generate all subvecs of the adapter and revcomp and insert them into the hashset
        for check_len in tolerance..(a.len() + 1) {
            let subadapter = Vec::from(&a[a.len() - check_len..]);
            adapters.insert(subadapter);

            let rev_subadapter = Vec::from(&rev[a.len() - check_len..]);
            adapters.insert(rev_subadapter);
        }
    }

    adapters
}

/// Parse the command line configuration from clap.
pub fn parse_config(args: &ArgMatches) -> MtsvResult<PrepConfig> {

    let num_threads = match args.value_of("NUM_THREADS") {
        Some(s) => s.parse::<usize>().expect("Invalid number entered for number of threads!"),
        None => unreachable!(),
    };

    let seg_len = match args.value_of("SEGMENT") {
        Some(l) => Some(l.parse::<usize>().expect("Invalid segment length provided")),
        None => None,
    };

    let (quality_min, quality_threshold) = match (args.value_of("QUALITY_MIN"),
                                                  args.value_of("QUALITY_THRESHOLD")) {
        (Some(m), Some(t)) => {
            let m = m.parse::<u8>().expect("Invalid quality minimum provided");
            let t = t.parse::<usize>().expect("Invalid quality # bases threshold provided");
            (Some(m), Some(t))
        },
        (None, None) => (None, None),
        (None, Some(_)) => panic!("Quality threshold was specified but not minimum quality."),
        (Some(_), None) => panic!("Quality minimum was specified but not # bases threshold."),
    };

    let (adapters, adapter_tolerance) = match (args.value_of("ADAPTER_FILE"),
                                               args.value_of("ADAPTER_TOLERANCE")) {
        (Some(f), Some(t)) => {
            let tolerance = t.parse::<usize>().expect("Invalid adapter tolerance provided");

            let adapters = try!(read_adapters(&Path::new(f), tolerance));

            (Some(adapters), Some(tolerance))
        },
        (None, None) => (None, None),
        (None, Some(_)) => panic!("Adapter tolerance provided, but not an adapter file."),
        (Some(_), None) => panic!("Adapter file provided, but not an adapter tolerance."),
    };

    let outfile = PathBuf::from(args.value_of("FASTA").unwrap());

    let mut infiles = Vec::new();
    info!("Parsing FASTQ files to determine minimum read length...");
    for p in args.values_of("FASTQ").unwrap() {
        info!("Parsing {}...", p);
        let p = PathBuf::from(p);
        let md = try!(read_fastq_metadata(&p));

        infiles.push((p, md));
    }

    let overall_min_read_len = infiles.iter().map(|&(_, ref md)| md.read_len).min().unwrap();

    let trim = if args.is_present("LCD") {
        TrimType::LcdFirstN(overall_min_read_len)
    } else if args.is_present("LCDQ") {
        TrimType::LcdQuality(overall_min_read_len)
    } else if args.is_present("SEGMENT") {
        TrimType::Segment(seg_len.unwrap())
    } else {
        unreachable!();
    };

    // validation:
    // make sure segment length is short enough to be satisfied by minimum read length
    // make sure quality threshold is short enough to be satisfied by minimum read length

    Ok(PrepConfig {
        trim: trim,
        min_quality: quality_min,
        quality_threshold: quality_threshold,
        adapters: adapters,
        adapter_tolerance: adapter_tolerance,
        num_threads: num_threads,
        infiles: infiles,
        outfile: outfile,
    })
}

/// Generate the clap Application for readprep.
pub fn prep_cli_app() -> App<'static, 'static> {
    App::new("readprep")
        .version(env!("CARGO_PKG_VERSION"))
        .author(env!("CARGO_PKG_AUTHORS"))
        .about("Read fragment quality control and homogenization tool (FASTQ -> FASTA).")
        .arg(Arg::with_name("NUM_THREADS")
            .short("t")
            .long("threads")
            .takes_value(true)
            .help("Number of worker threads to spawn.")
            .default_value("4"))
        .arg(Arg::with_name("VERBOSE")
            .short("v")
            .help("Include this flag to trigger debug-level logging."))
        .group(ArgGroup::with_name("TRIM")
            .arg("LCD")
            .arg("LCDQ")
            .arg("SEGMENT")
            .required(true))
        .arg(Arg::with_name("LCD")
            .long("lcd")
            .help("Enable LCD trim mode (takes first N bases of each read, where N = shortest \
                   read length in FASTQ files)."))
        .arg(Arg::with_name("LCDQ")
            .long("lcdqual")
            .help("Enable LCDQ trim mode (takes highest quality N bases of each read, where N = \
                   shortest read length in FASTQ files)."))
        .arg(Arg::with_name("SEGMENT")
            .long("segment")
            .help("Enable SEG trim mode (takes subsequent N length subsequences of each read).")
            .takes_value(true))
        .arg(Arg::with_name("ADAPTER_FILE")
            .long("adapters")
            .help("Path to file containing adapters, one per line.")
            .takes_value(true)
            .requires("ADAPTER_TOLERANCE"))
        .arg(Arg::with_name("ADAPTER_TOLERANCE")
            .long("adapter-tolerance")
            .help("Number of adapter characters to tolerate at start of reads.")
            .takes_value(true)
            .requires("ADAPTER_FILE"))
        .arg(Arg::with_name("QUALITY_MIN")
            .long("quality_min")
            .help("Minimum FASTQ quality to tolerate per base.")
            .takes_value(true)
            .requires("QUALITY_THRESHOLD"))
        .arg(Arg::with_name("QUALITY_THRESHOLD")
            .long("quality_threshold")
            .help("Maximum number of bases below minimum quality to tolerate per read.")
            .takes_value(true)
            .requires("QUALITY_MIN"))
        .arg(Arg::with_name("FASTA")
            .short("o")
            .long("out")
            .help("Path to desired output FASTA file.")
            .takes_value(true)
            .required(true))
        .arg(Arg::with_name("FASTQ")
            .help("Path(s) to FASTQ files to QC and collapse.")
            .takes_value(true)
            .multiple(true)
            .required(true)
            .display_order(1000)
            .validator(validate_fastq_file))
}

/// Check to make sure we can parse a FASTQ record out of a given path.
fn validate_fastq_file(s: String) -> Result<(), String> {
    let mut rdr = try!(fastq::Reader::from_file(Path::new(&s))
        .map_err(|e| format!("Unable to open {} for parsing: {:?}", &s, e)));
    let mut rec = fastq::Record::new();

    rdr.read(&mut rec).map_err(|e| format!("Unable to parse {} as a FASTQ file: {:?}", &s, e))
}

#[cfg(test)]
mod test {

    use std::collections::HashSet;
    use std::path::{Path, PathBuf};
    use super::*;

    #[test]
    fn segment() {
        let passed = "mtsv-readprep \
                      --segment 50 --out /dev/null \
                      tests/prep/sample1.fastq tests/prep/sample2.fastq";

        let expected = PrepConfig {
            trim: TrimType::Segment(50),
            num_threads: 4,
            min_quality: None,
            quality_threshold: None,
            adapter_tolerance: None,
            adapters: None,
            outfile: PathBuf::from("/dev/null"),
            infiles: vec![
                (PathBuf::from("tests/prep/sample1.fastq"),
                 FastqMetadata {
                    num_reads: 61927,
                    qual_encoding: QualityEncoding::Sanger,
                    read_len: 100,
                }),
                (PathBuf::from("tests/prep/sample2.fastq"),
                 FastqMetadata {
                    num_reads: 50402,
                    qual_encoding: QualityEncoding::Sanger,
                    read_len: 100,
                }),
            ],
        };

        let app = prep_cli_app();

        let args = app.get_matches_from_safe(passed.split(' ')).unwrap();
        let parsed = parse_config(&args).unwrap();

        assert_eq!(parsed, expected);
    }

    #[test]
    fn lcd() {
        let passed = "mtsv-readprep \
                      --lcd --out /dev/null \
                      tests/prep/sample1.fastq tests/prep/sample2.fastq";

        let expected = PrepConfig {
            trim: TrimType::LcdFirstN(100),
            num_threads: 4,
            min_quality: None,
            quality_threshold: None,
            adapter_tolerance: None,
            adapters: None,
            outfile: PathBuf::from("/dev/null"),
            infiles: vec![
                (PathBuf::from("tests/prep/sample1.fastq"),
                 FastqMetadata {
                    num_reads: 61927,
                    qual_encoding: QualityEncoding::Sanger,
                    read_len: 100,
                }),
                (PathBuf::from("tests/prep/sample2.fastq"),
                 FastqMetadata {
                    num_reads: 50402,
                    qual_encoding: QualityEncoding::Sanger,
                    read_len: 100,
                }),
            ],
        };

        let app = prep_cli_app();

        let args = app.get_matches_from_safe(passed.split(' ')).unwrap();
        let parsed = parse_config(&args).unwrap();

        assert_eq!(parsed, expected);
    }

    #[test]
    fn test_subadapters() {
        let expected = vec!["AAG",
                            "GAAG",
                            "GGAAG",
                            "CGGAAG",
                            "TCGGAAG",
                            "ATCGGAAG",
                            "GATCGGAAG",
                            "ATC",
                            "GATC",
                            "CGATC",
                            "CCGATC",
                            "TCCGATC",
                            "TTCCGATC",
                            "CTTCCGATC"];

        let expected = expected.into_iter().map(|a| Vec::from(a)).collect::<HashSet<Vec<u8>>>();
        let found = read_adapters(&Path::new("tests/prep/adapter.small"), 3).unwrap();

        assert_eq!(expected, found);
    }
}
