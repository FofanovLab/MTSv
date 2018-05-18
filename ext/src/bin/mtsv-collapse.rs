#[macro_use]
extern crate log;

extern crate clap;
extern crate mtsv;


use clap::{App, Arg};
use std::fs::File;
use std::io::{BufReader, BufWriter};

use mtsv::collapse::collapse_files;
use mtsv::util;

fn main() {
    let args = App::new("mtsv-collapse")
        .version(env!("CARGO_PKG_VERSION"))
        .author(env!("CARGO_PKG_AUTHORS"))
        .about("Tool for combining the output of multiple separate mtsv runs.")
        .arg(Arg::with_name("OUTPUT")
            .help("Path to write combined outupt file to.")
            .short("o")
            .long("output")
            .takes_value(true)
            .required(true))
        .arg(Arg::with_name("FILES")
            .index(1)
            .help("Path(s) to mtsv results files to collapse")
            .takes_value(true)
            .multiple(true)
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
    let files = args.values_of("FILES").unwrap().collect::<Vec<_>>();

    let mut infiles = Vec::new();

    // fail fast by open all the files to start
    info!("Opening output file...");
    let mut outfile = BufWriter::new(File::create(outpath).expect("Unable to create output file."));

    info!("Opening input files...");
    for f in files {
        let rdr = BufReader::new(File::open(f)
            .expect(&format!("Unable to open {} for reading.", f)));
        infiles.push(rdr);
    }

    match collapse_files(&mut infiles, &mut outfile) {
        Ok(()) => {
            info!("Successfully collapsed files. Output available in {}",
                  outpath)
        },
        Err(why) => panic!("Problem collapsing files: {}", why),
    }
}
