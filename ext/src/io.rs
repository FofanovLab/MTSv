//! Helper functions for serialization & deserialization.

use bincode::SizeLimit;
use bincode::rustc_serialize::{decode_from, encode_into};
use bio::io::fasta;
use error::*;
use index::{Database, TaxId};
use rustc_serialize::{Decodable, Encodable};
use std::collections::{BTreeMap, BTreeSet};
use std::fs::File;
use std::io;
use std::io::{BufRead, BufReader, BufWriter};
use std::path::Path;
use util::parse_read_header;

/// Parse an arbitrary `Decodable` type from a file path.
pub fn from_file<T>(p: &str) -> MtsvResult<T>
    where T: Decodable
{

    let f = try!(File::open(Path::new(p)));
    let mut reader = BufReader::new(f);
    Ok(try!(decode_from(&mut reader, SizeLimit::Infinite)))
}

/// Write an arbitrary `Encodable` type to a file path.
pub fn write_to_file<T>(t: &T, p: &str) -> MtsvResult<()>
    where T: Encodable
{

    let f = try!(File::create(Path::new(p)));
    let mut writer = BufWriter::new(f);
    Ok(try!(encode_into(t, &mut writer, SizeLimit::Infinite)))
}

/// Parse a FASTA database into a single map of all taxonomy IDs.
pub fn parse_fasta_db<R>(records: R) -> MtsvResult<Database>
    where R: Iterator<Item = io::Result<fasta::Record>>
{
    let mut taxon_map = BTreeMap::new();

    debug!("Parsing FASTA database file...");
    for record in records {
        let record = try!(record);

        let (gi, tax_id) = try!(parse_read_header(record.header()));

        let sequences = taxon_map.entry(tax_id).or_insert_with(|| vec![]);
        sequences.push((gi, record.seq().to_vec()));
    }

    Ok(taxon_map)
}

/// Return a lazy iterator which parses the findings of a mtsv-binner run.
///
/// The Option return type could indicate a few problems:
///
/// * There are an incorrect number of tokens after splitting on the colon separator
/// * One of the tax IDs isn't a valid unsigned integer
///
pub fn parse_findings<'a, R: BufRead + 'a>
    (s: R)
     -> Box<Iterator<Item = MtsvResult<(String, BTreeSet<TaxId>)>> + 'a> {
    // TODO: replace with -> impl Trait when stabilized

    // the BufRead::lines function handles lazily splitting on lines for us
    Box::new(s.lines().map(|l| {
        l.map_err(|e| MtsvError::from(e)).and_then(|l| {
            let l = l.trim();
            // split from the right in case someone put colons in the read ID
            let mut halves = l.rsplitn(2, ':');

            let mut hits = BTreeSet::new();

            // the first split iteration will always return something, even if it's empty
            let taxids = halves.next().unwrap().split(',');

            // parse each taxid (comma separated), returning None if it fails
            for taxid_raw in taxids {
                let taxid = match taxid_raw.parse::<TaxId>() {
                    Ok(id) => id,
                    Err(_) => return Err(MtsvError::InvalidInteger(taxid_raw.to_string())),
                };

                hits.insert(taxid);
            }

            // since we're parsing from the right of each line, the read ID is the second token
            let read_id = match halves.next() {
                Some(r) => {
                    if r.len() > 0 {
                        r.to_string()
                    } else {
                        return Err(MtsvError::InvalidHeader(l.to_string()));
                    }
                },
                None => return Err(MtsvError::InvalidHeader(l.to_string())),
            };

            Ok((read_id, hits))
        })
    }))
}

#[cfg(test)]
mod test {

    use ::binner::write_single_line;
    use ::index::TaxId;

    use mktemp::Temp;

    use rand::{Rng, XorShiftRng};
    use std::collections::{BTreeMap, BTreeSet};
    use std::io::{BufReader, Cursor};
    use std::iter::FromIterator;
    use super::*;

    fn roundtrip(findings: Vec<(String, BTreeSet<TaxId>)>) {

        let mut buf = Vec::new();

        for &(ref header, ref matches) in &findings {
            write_single_line(header, &matches, &mut buf).unwrap();
        }

        let results = parse_findings(Cursor::new(buf));

        let mut expected = findings.into_iter();

        for res in results {
            let (found_head, found_matches) = res.unwrap();
            let (expected_head, expected_matches) = expected.next().unwrap();
            assert_eq!(found_head, expected_head);
            assert_eq!(found_matches, expected_matches);
        }
    }

    #[test]
    fn roundtrip_single() {
        let header = String::from("raldkjfasdlkfj");
        let mut matches = BTreeSet::new();
        matches.insert(TaxId(2093874));
        matches.insert(TaxId(12334));
        matches.insert(TaxId(65198));
        matches.insert(TaxId(1309579821));
        matches.insert(TaxId(241324));

        roundtrip(vec![(header, matches)]);
    }

    #[test]
    fn roundtrip_many() {
        let mut rng = XorShiftRng::new_unseeded();

        let num_findings: usize = rng.gen_range(500, 1_000);

        let mut findings = Vec::with_capacity(num_findings);

        for _ in 0..num_findings {
            let header_len: usize = rng.gen_range(1, 100);
            let num_matches: usize = rng.gen_range(1, 1_000);

            let header: String = rng.gen_ascii_chars().take(header_len).collect();
            let mut matches = BTreeSet::new();

            for _ in 0..num_matches {
                matches.insert(TaxId(rng.gen()));
            }

            findings.push((header, matches));
        }

        roundtrip(findings);
    }

    #[test]
    fn parsing_positive() {
        let working = String::from("r1234:1,2,3
r12345:5,7,3
asldkfj:3,4,5,6")
            .into_bytes();

        let expected = {
            let mut e = BTreeMap::new();
            e.insert(String::from("r1234"),
                     BTreeSet::from_iter(vec![TaxId(1), TaxId(2), TaxId(3)].into_iter()));

            e.insert(String::from("r12345"),
                     BTreeSet::from_iter(vec![TaxId(5), TaxId(7), TaxId(3)].into_iter()));

            e.insert(String::from("asldkfj"),
                     BTreeSet::from_iter(vec![TaxId(3), TaxId(4), TaxId(5), TaxId(6)].into_iter()));

            e
        };

        let mut results = BTreeMap::new();

        for res in parse_findings(working.as_slice()) {
            let (read_header, hits) = res.unwrap();
            results.insert(read_header, hits);
        }

        assert_eq!(expected, results);
    }

    #[test]
    #[should_panic]
    fn missing_ids() {
        let bad = String::from(":");
        let bad = BufReader::new(Cursor::new(bad.as_bytes()));

        for i in parse_findings(bad) {
            i.unwrap();
        }
    }

    #[test]
    #[should_panic]
    fn invalid_ids() {
        let bad = String::from("r12345:abc,def,ghi");
        let bad = BufReader::new(Cursor::new(bad.as_bytes()));

        for i in parse_findings(bad) {
            i.unwrap();
        }
    }

    #[test]
    #[should_panic]
    fn no_read_header() {
        let bad = String::from("123,456,789");
        let bad = BufReader::new(Cursor::new(bad.as_bytes()));

        for i in parse_findings(bad) {
            i.unwrap();
        }
    }

    quickcheck! {
        fn io_helpers(map: BTreeMap<String, String>) -> bool {
            let outfile = Temp::new_file().unwrap();
            let outfile = outfile.to_path_buf();
            let outfile = outfile.to_str().unwrap();

            write_to_file(&map, outfile).unwrap();
            let from_file = from_file(outfile).unwrap();

            map == from_file
        }
    }
}
