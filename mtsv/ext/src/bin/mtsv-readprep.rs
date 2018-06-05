#[macro_use]
extern crate log;

extern crate bio;
extern crate clap;
extern crate rustc_serialize;
extern crate mtsv;

use mtsv::prep::*;
use mtsv::prep_config::*;
use mtsv::util;

fn main() {
    let args = prep_cli_app().get_matches();

    // setup logger
    util::init_logging(if args.is_present("VERBOSE") {
        log::LogLevelFilter::Debug
    } else {
        log::LogLevelFilter::Info
    });

    let config = match parse_config(&args) {
        Ok(c) => c,
        Err(why) => panic!("Unable to generate config: {:?}", why),
    };

    println!("{:?}", &config);

    let exit_code = match run_prep(&config) {
        Ok(()) => {
            info!("Prep finished!");
            0
        },
        Err(why) => {
            error!("Prep unsuccessful! ({:?})", why);
            1
        },
    };

    std::process::exit(exit_code);
}
