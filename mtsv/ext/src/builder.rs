//! Build metagenomic index for binning queries.

use bio::io::fasta;

use error::*;
use index::MGIndex;
use io::{parse_fasta_db, write_to_file};
use std::io;

/// Build and write the metagenomic index to disk.
///
/// The actual construction logic is in `mtsv::index::MGIndex`, this just handles the I/O and
/// parsing.
pub fn build_and_write_index<R>(records: R,
                                index_path: &str,
                                sample_interval: u32)
                                -> MtsvResult<()>
    where R: Iterator<Item = io::Result<fasta::Record>>
{
    let taxon_map = try!(parse_fasta_db(records));

    info!("File parsed, building index...");
    let index = MGIndex::new(taxon_map, sample_interval);

    info!("Writing index to file...");
    try!(write_to_file(&index, index_path));

    Ok(())
}

#[cfg(test)]
mod test {
    use bio::io::fasta::Reader;
    use mktemp::Temp;
    use std::io::Cursor;
    use super::build_and_write_index;

    #[test]
    fn success() {
        let reference = ">123-456
TGTCTTAATGATAAAAATTGTTACAAACAGTTTAACATATTTAGCTACCTATTTTGCATATAAAAAACATGCTTGCATACACTATGCAATAAAAATTACAAATTTATATATGATACCACTATGCTTGCTTATCTCTATAGCGCCATTGATACACATTTTTAAATATCTATACTGCCGTTAGAATTTTATCATGTCTTAATTTTCATTAAATATTAATTACTTCATTTTATATAAACCAACAAAAACCCCCTCACTACTATGCAAGTGAGAGGTTATGTTGATGTGCTTTATTTTCAT
\
                         >124-456
TTTCACCTAGTACATTAAATACACGACCTAATGTTTCGTCACCAACAGGTACACTAATTTCTTTGCCTGTATCTTTTACATCCATGCCTCTTTGGACACCATCAGTTGAATCCATCGCAATTGTACGAACAACGTCGTCACCTAATTGCAGCGCAACTTCTAATGTTAGTTGTATTGTACCTTCTTCTTTAGGCACATCAATAACCAAGGCGTTATTAATTTTAGGAACTTCGTTATGTTCAAATCGAACATCAATTACAGGACCCATAACTTGAGTTACACGGCCAATTCCCATGCTATTTTCCTCCTTTAAATATTATTCAAGCGCTGCGGAACCACCAACAATTTCAGTAATTTGTTGCGTAATTTCTGCTTGTCTCGCTCTGTTATATTCTA
\
                         >908-678
AAAACACATATTTTCAAATCTAGTAAATATTAAATCTACTCTTGACGATTGCACCAATGCTACGCGATATAGATATCCACTAAAAACATACGTAATCATAACCATCATTGTTAGAAACAAAATTATTTCCATGATAACCCTCACTTAATATATTTCTAAAATTTTTCACTACGAATTAAGGCATAAAATAAATACAAAACTAATGCAATAACTACCAGTAATAAAACGATGAGCATTGCCATAACC";

        let records = Reader::new(Cursor::new(reference.as_bytes())).records();
        let outfile = Temp::new_file().unwrap();
        let outfile_path = outfile.to_path_buf();
        let outfile_str = outfile_path.to_str().unwrap();

        build_and_write_index(records, outfile_str, 32).unwrap();

        assert!(outfile_path.exists());
        assert!(outfile_path.is_file());

        let metadata = outfile_path.metadata().unwrap();

        assert!(metadata.len() > reference.len() as u64);
    }

    #[test]
    #[should_panic]
    fn fail_empty_header() {
        let reference = ">
TGTCTTAATGATAAAAATTGTTACAAACAGTTTAACATATTTAGCTACCTATTTTGCATATAAAAAACATGCTTGCATACACTATGCAATAAAAATTACAAATTTATATATGATACCACTATGCTTGCTTATCTCTATAGCGCCATTGATACACATTTTTAAATATCTATACTGCCGTTAGAATTTTATCATGTCTTA
\
                         >124-456
TTTCACCTAGTACATTAAATACACGACCTAATGTTTCGTCACCAACAGGTACACTAATTTCTTTGCCTGTATCTTTTACATCCATGCCTCTTTGGACACCATCAGTTGAATCCATCGCAATTGTACGAACAACGTCGTCACCTAATTGCAGCGCAACTTCTAATGTTAGTTGTATTGTACCTTCTTCTTTAGGCACATCAATAACCAAGGCGTTATTAATTTTAGGAACTTCGTTATGTTCAAATCGAACATCAATTACAGGACCCATAACTTGAGTTACACGGCCAATTCCCATGC";

        let records = Reader::new(Cursor::new(reference.as_bytes())).records();
        let outfile = Temp::new_file().unwrap();
        let outfile_path = outfile.to_path_buf();
        let outfile_str = outfile_path.to_str().unwrap();

        build_and_write_index(records, outfile_str, 32).unwrap();
    }
}
