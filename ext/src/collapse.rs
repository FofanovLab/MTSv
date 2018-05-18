//! Collapse multiple mtsv results/findings files into a single one.

use binner::write_single_line;
use error::*;
use io::parse_findings;
use std::collections::{BTreeMap, BTreeSet};
use std::io::{BufRead, Write};

/// Given a list of mtsv results file paths, collapse into a single one.
pub fn collapse_files<R, W>(files: &mut [R], write_to: &mut W) -> MtsvResult<()>
    where R: BufRead,
          W: Write
{
    let mut results = BTreeMap::new();

    for ref mut r in files {

        for res in parse_findings(r) {
            let (readid, hits) = try!(res);

            results.entry(readid).or_insert(BTreeSet::new()).extend(hits.into_iter());
        }
    }

    info!("All input files parsed and collapsed, writing to disk...");
    for (header, hits) in results.iter() {
        try!(write_single_line(header, hits, write_to));
    }

    Ok(())
}

#[cfg(test)]
mod test {
    use std::io::Cursor;
    use super::*;

    #[test]
    fn simple_collapse() {
        let a = "a:1,2,3,4,5
b:1,2,3,4";
        let b = "b:3,4,5,6,7
a:8,9,10,100";
        let c = "c:2,3,4,5";

        let mut buf = Vec::new();
        let mut buf2 = Vec::new();

        let mut infiles = vec![Cursor::new(a), Cursor::new(b), Cursor::new(c)];
        let mut infiles2 = vec![Cursor::new(b), Cursor::new(c), Cursor::new(a)];

        collapse_files(&mut infiles, &mut buf).unwrap();
        collapse_files(&mut infiles2, &mut buf2).unwrap();

        let buf_str = String::from_utf8(buf).unwrap();
        let buf2_str = String::from_utf8(buf2).unwrap();

        assert_eq!(buf_str, buf2_str);

        let expected = "a:1,2,3,4,5,8,9,10,100
b:1,2,3,4,5,6,7
c:2,3,4,5
";

        assert_eq!(expected, &buf_str);
    }
}
