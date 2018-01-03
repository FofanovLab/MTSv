#[macro_use]
extern crate log;

extern crate bio;
extern crate clap;
extern crate flate2;
extern crate tar;
extern crate mtsv;


use clap::{App, Arg};
use std::fs::File;
use std::io::{BufReader, BufWriter};

use mtsv::io::{from_file, parse_findings};
use mtsv::tax_tree::{LcaSetting, LogicalSetting, TreeWithIndices};
use mtsv::util;

fn main() {

    let args = App::new("mtsv-inform")
        .version(env!("CARGO_PKG_VERSION"))
        .author(env!("CARGO_PKG_AUTHORS"))
        .about("Postprocessor for mtsv results to determine which reads are \"informative.\"")
        .arg(Arg::with_name("INPUT")
            .short("i")
            .long("input")
            .help("Path to mtsv results file.")
            .takes_value(true)
            .required(true))
        .arg(Arg::with_name("OUTPUT")
            .short("o")
            .long("output")
            .help("Output path to write \"informativeness\" results.")
            .takes_value(true)
            .required(true))
        .arg(Arg::with_name("INDEX")
            .short("x")
            .long("index")
            .help("Path to index built by mtsv-tree-build utility.")
            .takes_value(true)
            .required(true))
        .arg(Arg::with_name("GENUS")
            .long("genus")
            .help("Enable to search for a common GENUS among hits for a read (takes priority \
                   over LCA search when a family exists for a taxonomic ID).")
            .conflicts_with("FAMILY"))
        .arg(Arg::with_name("FAMILY")
            .long("family")
            .help("Enable to search for a common FAMILY among hits for a read (takes priority \
                   over LCA search when a family exists for a taxonomic ID).")
            .conflicts_with("GENUS"))
        .arg(Arg::with_name("LCA")
            .long("lca")
            .help("Height at which the search will attempt to find a common ancestor among the \
                   hits for a read.")
            .possible_values(&["0", "1", "2", "3"])
            .required(true)
            .takes_value(true))
        .arg(Arg::with_name("NUM_THREADS")
            .short("t")
            .long("threads")
            .takes_value(true)
            .help("Number of worker threads to spawn.")
            .default_value("4"))
        .arg(Arg::with_name("VERBOSE")
            .short("v")
            .help("Include this flag to trigger debug-level logging."))
        .get_matches();

    let num_threads = match args.value_of("NUM_THREADS") {
        Some(s) => s.parse::<usize>().expect("Invalid number entered for number of threads!"),
        None => unreachable!(),
    };

    let lca = match args.value_of("LCA").unwrap().parse::<u8>().unwrap() {
        0 => LcaSetting::Zero,
        1 => LcaSetting::One,
        2 => LcaSetting::Two,
        _ => LcaSetting::Three,
    };

    let logical = if args.is_present("GENUS") {
        LogicalSetting::Genus
    } else if args.is_present("FAMILY") {
        LogicalSetting::Family
    } else {
        LogicalSetting::Na
    };

    // setup logger
    util::init_logging(if args.is_present("VERBOSE") {
        log::LogLevelFilter::Debug
    } else {
        log::LogLevelFilter::Info
    });

    let input_path = args.value_of("INPUT").unwrap();
    let index_path = args.value_of("INDEX").unwrap();
    let output_path = args.value_of("OUTPUT").unwrap();

    info!("Opening files...");
    let input_rdr = match File::open(input_path) {
        Ok(f) => BufReader::new(f),
        Err(why) => panic!("Unable to open input file: {}", why),
    };

    let mut output_wtr = match File::create(output_path) {
        Ok(f) => BufWriter::new(f),
        Err(why) => panic!("Unable to create output file: {}", why),
    };

    info!("Deserializing taxonomic index...");
    let index = match from_file::<TreeWithIndices>(index_path) {
        Ok(i) => i,
        Err(why) => panic!("Unable to deserialize taxonomic index: {}", why),
    };

    info!("Finding informative reads...");
    index.find_and_write_informatives(parse_findings(input_rdr),
                                      num_threads,
                                      lca,
                                      logical,
                                      &mut output_wtr);

    info!("Successfully written to disk, exiting.");
}
