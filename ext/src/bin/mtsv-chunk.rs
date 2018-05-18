//! Split a FASTA reference database file into chunks well-suited for vedro's index generation.

#[macro_use]
extern crate log;

extern crate bio;
extern crate clap;
extern crate mtsv;

use bio::io::fasta;
use clap::{App, Arg};
use std::path::Path;
use mtsv::chunk::write_db_chunks;
use mtsv::io::parse_fasta_db;
use mtsv::util;

fn main() {
    let args = App::new("vedro-chunk")
        .version(env!("CARGO_PKG_VERSION"))
        .author(env!("CARGO_PKG_AUTHORS"))
        .about("Split a FASTA reference database into chunks for index generation.")
        .arg(Arg::with_name("OUTPUT")
            .help("Folder path to write split outupt files to.")
            .short("o")
            .long("output")
            .takes_value(true)
            .required(true))
        .arg(Arg::with_name("INPUT")
            .help("Path(s) to vedro results files to collapse")
            .short("i")
            .long("input")
            .takes_value(true)
            .required(true))
        .arg(Arg::with_name("SIZE_GB")
            .help("Chunk size (in gigabytes).")
            .short("g")
            .long("gb")
            .default_value("1.0")
            .takes_value(true)
            .required(true))
        .arg(Arg::with_name("VERBOSE")
            .short("v")
            .help("Include this flag to trigger debug-level logging."))
        .get_matches();


    // setup logger
    util::init_logging(if args.is_present("VERBOSE") {
        log::LogLevelFilter::Debug
    } else {
        log::LogLevelFilter::Info
    });

    let outpath = args.value_of("OUTPUT").unwrap();
    let database = args.value_of("INPUT").unwrap();

    let base_name = match Path::new(database).file_stem().and_then(|s| s.to_str()) {
        Some(f) => f,
        None => panic!("{} is not a valid file.", database),
    };

    info!("Will write files as {}_0.fasta, {}_1.fasta, etc. in output path.",
          base_name,
          base_name);

    let chunks_gb = args.value_of("SIZE_GB").unwrap();
    let chunks_gb = chunks_gb.parse::<f32>()
        .unwrap_or_else(|_| {
            chunks_gb.parse::<u32>()
                .map(|u| u as f32)
                .expect("Unable to parse SIZE_GB as a valid number.")
        });

    if chunks_gb <= 0.0 {
        panic!("Unable to write negatively sized database chunks.");
    }

    let records = match fasta::Reader::from_file(database) {
        Ok(r) => r.records(),
        Err(why) => {
            panic!("Unable to open database file: {:?}", why);
        },
    };

    info!("Parsing database file...");

    let database = match parse_fasta_db(records) {
        Ok(d) => d,
        Err(why) => panic!("Unable to parse FASTA DB file: {:?}", why),
    };

    info!("Done parsing database file.");

    match write_db_chunks(&database, &base_name, Path::new(outpath), chunks_gb) {
        Ok(paths) => info!("Finished writing chunks to: {:#?}", paths),
        Err(why) => panic!("Unable to write DB chunks: {:?}", why),
    }
}
