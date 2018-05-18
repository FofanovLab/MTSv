//! Utilities for chunking database files.

use error::*;
use index::Database;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

/// Write database sequences to a series of files
pub fn write_db_chunks(records: &Database,
                       base_filename: &str,
                       out_path: &Path,
                       chunk_gbs: f32)
                       -> MtsvResult<Vec<PathBuf>> {

    if !out_path.is_dir() {
        return Err(MtsvError::MissingFile(format!("{} is not a directory",
                                                   out_path.to_string_lossy())));
    }

    let mut chunk_num = 0;
    let mut bytes_written = 0;
    let target_size = (chunk_gbs * 1_000_000_000.0) as usize;

    let mut written_paths = Vec::new();
    let mut chunk_path = out_path.to_path_buf();

    // add the actual filename
    chunk_path.push(&format!("{}_{}.fasta", base_filename, chunk_num));
    written_paths.push(chunk_path.clone());

    let mut writer = BufWriter::new(try!(File::create(&chunk_path)));

    info!("Writing to {:?}...", &chunk_path);

    for (tax_id, seqs) in records {
        let tid_str = tax_id.0.to_string();

        for &(gi, ref sequence) in seqs {
            let gi_str = gi.0.to_string();

            bytes_written += try!(writer.write(b">"));
            bytes_written += try!(writer.write(gi_str.as_bytes()));
            bytes_written += try!(writer.write(b"-"));
            bytes_written += try!(writer.write(tid_str.as_bytes()));
            bytes_written += try!(writer.write(b"\n"));
            bytes_written += try!(writer.write(&sequence));
            bytes_written += try!(writer.write(b"\n"));

            if bytes_written >= target_size {
                // we need to set up a new writer before we write things
                bytes_written = 0;
                chunk_num += 1;
                chunk_path.pop();
                chunk_path.push(&format!("{}_{}.fasta", base_filename, chunk_num));
                written_paths.push(chunk_path.clone());
                writer = BufWriter::new(try!(File::create(&chunk_path)));

                info!("Writing to {:?}...", &chunk_path);
            }
        }
    }

    Ok(written_paths)
}

#[cfg(test)]
mod test {
    use bio::io::fasta;
    use index::{Database, random_database};
    use io::parse_fasta_db;
    use mktemp::Temp;
    use std::fmt::Debug;
    use std::path::Path;
    use super::*;

    fn collect_chunks<P: AsRef<Path> + Debug>(paths: &[P]) -> Database {

        let mut overall = Database::new();

        for path in paths {
            println!("reading {:?}", path);
            let records = fasta::Reader::from_file(path).unwrap().records();
            let database = parse_fasta_db(records).unwrap();

            for (tax_id, seqs) in database {
                overall.entry(tax_id).or_insert_with(Vec::new).extend(seqs);
            }
        }

        overall
    }

    #[test]
    fn chunk_roundtrip() {
        let db = random_database(100, 200, 500, 10_000);

        let dir = Temp::new_dir().unwrap();
        let dir = dir.to_path_buf();

        let chunks = write_db_chunks(&db, "tmp_fasta", &dir, 0.001).unwrap();

        let expected = collect_chunks(&chunks);

        assert_eq!(db, expected);
    }
}
