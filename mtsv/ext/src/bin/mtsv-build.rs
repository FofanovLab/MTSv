#[macro_use]
extern crate log;

extern crate bio;
extern crate clap;
extern crate mtsv;


use bio::io::fasta;
use clap::{App, Arg};
use std::path::Path;
use mtsv::builder;
use mtsv::util;

fn main() {

    let args = App::new("mtsv-build")
        .version(env!("CARGO_PKG_VERSION"))
        .author(env!("CARGO_PKG_AUTHORS"))
        .about("Index construction for mtsv metagenomics binning tool.")
        .arg(Arg::with_name("FASTA")
            .short("f")
            .long("fasta")
            .help("Path to FASTA database file.")
            .takes_value(true)
            .required(true))
        .arg(Arg::with_name("INDEX")
            .short("i")
            .long("index")
            .help("Absolute path to mtsv index file.")
            .takes_value(true)
            .required(true))
        .arg(Arg::with_name("VERBOSE")
            .short("v")
            .help("Include this flag to trigger debug-level logging."))
        .arg(Arg::with_name("FM_SAMPLE_INTERVAL")
            .long("sample-interval")
            .takes_value(true)
            .help("Sampling interval for index generation. Smaller = more memory usage, \
                   slightly  faster queries. Larger = less memory usage slightly slower queries.")
            .default_value("32"))
        .get_matches();


    // setup logger
    util::init_logging(if args.is_present("VERBOSE") {
        log::LogLevelFilter::Debug
    } else {
        log::LogLevelFilter::Info
    });

    let fasta_path = args.value_of("FASTA").unwrap();
    let index_path = args.value_of("INDEX").unwrap();

    let exit_code = {

        let fm_index_interval = match args.value_of("FM_SAMPLE_INTERVAL") {
            Some(s) => s.parse::<u32>().expect("Invalid index sample interval entered!"),
            None => unreachable!(),
        };

        debug!("Opening FASTA database file...");
        let records = fasta::Reader::from_file(Path::new(fasta_path))
            .expect("Unable to open FASTA database for parsing.")
            .records();

        match builder::build_and_write_index(records, index_path, fm_index_interval) {
            Ok(_) => {
                info!("Done building and writing index!");
                0
            },
            Err(why) => {
                error!("Error building index: {}", why);
                1
            },
        }
    };

    std::process::exit(exit_code);
}
