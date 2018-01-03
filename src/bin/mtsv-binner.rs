#[macro_use]
extern crate log;

extern crate clap;

extern crate mtsv;

use clap::{App, Arg};

use mtsv::binner;
use mtsv::util;

fn main() {

    let args = App::new("mtsv")
        .version(env!("CARGO_PKG_VERSION"))
        .author(env!("CARGO_PKG_AUTHORS"))
        .about("Metagenomics binning tool.")
        .arg(Arg::with_name("FASTA")
            .short("f")
            .long("fasta")
            .help("Absolute path to FASTA query file.")
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
        .arg(Arg::with_name("RESULTS_PATH")
            .short("r")
            .long("results")
            .takes_value(true)
            .help("Path to write results file to."))
        .arg(Arg::with_name("NUM_THREADS")
            .short("t")
            .long("threads")
            .takes_value(true)
            .help("Number of worker threads to spawn.")
            .default_value("4"))
        .arg(Arg::with_name("EDIT_TOLERANCE")
            .short("e")
            .long("edits")
            .takes_value(true)
            .help("Edit distance to tolerate in matched reference sites")
            .default_value("3"))
        .arg(Arg::with_name("SEED_SIZE")
            .long("seed-size")
            .takes_value(true)
            .help("Set to override inital exact match query size.")
            .default_value("20"))
        .arg(Arg::with_name("SEED_GAP")
            .long("seed-gap")
            .takes_value(true)
            .help("Set to override gap between seeds used for initial exact match.")
            .default_value("3"))
        .arg(Arg::with_name("MIN_SEEDS")
            .long("min-seeds")
            .takes_value(true)
            .help("Set to override minimum number of seeds to perform alignment of a candidate \
                   site.")
            .default_value("2"))
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
        let results_path = args.value_of("RESULTS_PATH");

        let num_threads = match args.value_of("NUM_THREADS") {
            Some(s) => s.parse::<usize>().expect("Invalid number entered for number of threads!"),
            None => unreachable!(),
        };

        let edit_tolerance = match args.value_of("EDIT_TOLERANCE") {
            Some(s) => s.parse::<u32>().expect("Invalid edit distance entered!"),
            None => unreachable!(),
        };

        let seed_size = match args.value_of("SEED_SIZE") {
            Some(s) => {
                let seed_size = s.parse::<usize>().expect("Invalid seed size entered!");

                if seed_size < 16 {
                    warn!("Seed size may be small enough that it causes performance issues.");
                } else if seed_size > 24 {
                    warn!("Seed size may be large enough that significant results are ignored.");
                }

                seed_size
            },
            None => panic!("Missing parameter: seed_size"),
        };

        let seed_gap = match args.value_of("SEED_GAP") {
            Some(s) => {
                let seed_gap = s.parse::<usize>().expect("Invalid seed gap entered!");

                if seed_gap < 3 {
                    warn!("Seed gap may be small enough that it causes performance issues.");
                } else if seed_gap > 10 {
                    warn!("Seed gap may be large enough that significant results are ignored.");
                }

                seed_gap
            },
            None => panic!("Missing parameter: seed_gap"),
        };

        let min_seeds = match args.value_of("MIN_SEEDS") {
            Some(s) => {
                let min_seeds = s.parse::<usize>().expect("Invalid min. # of seeds entered!");

                if min_seeds < 2 {
                    warn!("Performance may be significantly slowed by aligning candidate sites \
                           with that few seeds found.");
                } else if min_seeds > 4 {
                    warn!("Min. # of seeds may be high enough that significant results are \
                           ignored.");
                }

                min_seeds
            },
            None => panic!("Missing parameter: seed_gap"),
        };

        if results_path.is_none() {
            error!("No results path provided!");
            3
        } else {
            let results_path = results_path.unwrap();
            match binner::get_and_write_matching_bin_ids(fasta_path,
                                                         index_path,
                                                         results_path,
                                                         num_threads,
                                                         edit_tolerance,
                                                         seed_size,
                                                         seed_gap,
                                                         min_seeds) {
                Ok(_) => 0,
                Err(why) => {
                    error!("Error running query: {}", why);
                    2
                },
            }
        }
    };

    std::process::exit(exit_code);
}
