#[macro_use]
extern crate log;

extern crate bio;
extern crate clap;
extern crate flate2;
extern crate tar;
extern crate mtsv;


use clap::{App, Arg};
use flate2::bufread::GzDecoder;
use std::fs::File;
use std::io;
use std::io::{Cursor, Read};
use tar::Archive;

use mtsv::error::*;
use mtsv::io::write_to_file;
use mtsv::tax_tree::TreeWithIndices;
use mtsv::util;

fn main() {

    let args = App::new("mtsv-tree-build")
        .version(env!("CARGO_PKG_VERSION"))
        .author(env!("CARGO_PKG_AUTHORS"))
        .about("Index construction for mtsv-inform read informativeness tool.")
        .arg(Arg::with_name("NCBI_TREE")
            .short("d")
            .long("dump")
            .help("Path to NCBI taxdump.tar.gz file which matches the given FASTA file.")
            .takes_value(true)
            .required(true))
        .arg(Arg::with_name("INDEX")
            .short("i")
            .long("index")
            .help("Output path to mtsv-inform index file.")
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

    let dump_path = args.value_of("NCBI_TREE").unwrap();
    let index_path = args.value_of("INDEX").unwrap();

    info!("Retrieving node dump from archive...");
    let dump = match get_node_dump_from_tar(&dump_path) {
        Ok(d) => Cursor::new(d),
        Err(why) => panic!("Problem retrieving node dump from archive ({:?})", why),
    };

    info!("Parsing NCBI dump for taxonomic tree...");
    match TreeWithIndices::from_node_dump(dump) {
        Ok(tree) => {
            info!("Parsing successful, writing tree to disk...");

            match write_to_file(&tree, &index_path) {
                Ok(()) => info!("Successfully wrote taxonomy tree index to disk"),
                Err(why) => panic!("Problem writing to disk: {}", why),
            }
        },
        Err(why) => {
            panic!("{:?}", why);
        },
    }
}

fn get_node_dump_from_tar(p: &str) -> MtsvResult<String> {
    let mut archive = Archive::new(try!(GzDecoder::new(io::BufReader::new(try!(File::open(p))))));

    for entry in try!(archive.entries()) {
        let mut entry = try!(entry);

        let mut bytes = Vec::new();
        try!(entry.read_to_end(&mut bytes));

        let path = try!(entry.path());
        let path_str = path.to_string_lossy();

        match &*path_str {
            "nodes.dmp" => {
                info!("Found nodes.dmp, lexing tax IDs...");
                return Ok(String::from_utf8_lossy(&bytes).into_owned());
            },
            _ => debug!("unneeded file: {}", path_str),
        }
    }

    Err(MtsvError::MissingFile("nodes.dmp".to_string()))
}
