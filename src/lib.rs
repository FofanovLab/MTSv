//! `vedro` is an umbrella project for several metagenomic binning tools.
//!
//! Currently:
//!
//! * `vedro-readprep`: a quality-control and deduplication tool which converts multiple FASTQ
//! files to a single FASTA file
//! * `vedro-binner`: a metagenomic binning tool that is sensitive to SNPs and other edits
//! * `vedro-build`: an index-construction tool to provide fast-lookup data structures for the
//! binner
//! * `vedro-collapse`: a simple tool to combine the output of several runs of the binner
//! * `vedro-tree-build`: an index-construction tool for NCBI taxonomic tree dumps -- these are
//! used in deciding which results are "informative"
//! * `vedro-inform`: takes findings from the binner/collapser, and decides whether they point to
//! the presence of a particular species/family/genus based on a user-specified level of precision
//!
//! All of these CLI tools depend on functionality in `vedro` -- they are defined as CLI-param
//! parsers which then call into library functionality.

#![warn(missing_docs)]

#[macro_use]
extern crate log;

extern crate bincode;
extern crate bio;
extern crate chrono;
extern crate clap;
extern crate cue;
extern crate daggy;
extern crate env_logger;
extern crate itertools;
extern crate rustc_serialize;
extern crate ssw;
extern crate stopwatch;

#[cfg(test)]
extern crate mktemp;

#[cfg(test)]
#[macro_use]
extern crate quickcheck;

#[cfg(test)]
extern crate rand;

pub mod align;
pub mod binner;
pub mod builder;
pub mod chunk;
pub mod collapse;
pub mod error;
pub mod index;
pub mod io;
pub mod prep;
pub mod prep_config;
pub mod tax_tree;
pub mod util;
