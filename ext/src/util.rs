//! Miscellaneous.

use chrono::Local;
use env_logger::LogBuilder;
use error::*;
use index::{Gi, TaxId};
use log::{LogLevelFilter, LogRecord};

/// Initialize the program-wide logger to write to stdout with timestamps.
pub fn init_logging(level: LogLevelFilter) {
    let mut builder = LogBuilder::new();

    builder.filter(None, level)
        .format(|record: &LogRecord| {
            format!("[{} {} {}] {}",
                    record.level(),
                    Local::now().format("%Y-%m-%d %H:%M:%S%.3f"),
                    record.location().module_path(),
                    record.args())
        });

    let _ = builder.init();
}

/// Parse a reference sequence's read header in the format expected by mtsv: `ACCESSION-TAXID`.
pub fn parse_read_header(h: &str) -> MtsvResult<(Gi, TaxId)> {
    let mut tokens = h.split('-');

    let gi = match tokens.next() {
        Some(t) => {
            match t.parse::<Gi>() {
                Ok(t) => t,
                Err(_) => return Err(MtsvError::InvalidInteger(t.to_owned())),
            }
        },
        None => return Err(MtsvError::InvalidHeader(String::from(h))),
    };

    let tax_id = match tokens.next() {
        Some(t) => {
            match t.parse::<TaxId>() {
                Ok(t) => t,
                Err(_) => return Err(MtsvError::InvalidInteger(t.to_owned())),
            }
        },
        None => return Err(MtsvError::InvalidHeader(String::from(h))),
    };

    if let None = tokens.next() {
        Ok((gi, tax_id))
    } else {
        // there's a second dash -- not the format we're expecting
        Err(MtsvError::InvalidHeader(String::from(h)))
    }
}

#[cfg(test)]
mod test {
    use index::{Gi, TaxId};

    use log::LogLevelFilter;
    use super::{init_logging, parse_read_header};

    #[test]
    fn lines_for_the_line_throne() {
        init_logging(LogLevelFilter::Debug);
    }

    #[test]
    fn success() {
        let (found_gi, found_tax) = parse_read_header("12345-908").unwrap();

        assert_eq!(found_gi, Gi(12345));
        assert_eq!(found_tax, TaxId(908));
    }

    #[test]
    #[should_panic]
    fn fail_empty_nodash() {
        let _ = parse_read_header("").unwrap();
    }

    #[test]
    #[should_panic]
    fn fail_empty() {
        let _ = parse_read_header("-").unwrap();
    }

    #[test]
    #[should_panic]
    fn fail_decimal_gi() {
        let _ = parse_read_header("1.0-543").unwrap();
    }

    #[test]
    #[should_panic]
    fn fail_decimal_taxid() {
        let _ = parse_read_header("654981-1.071").unwrap();
    }

    #[test]
    #[should_panic]
    fn fail_extra() {
        let _ = parse_read_header("1-2-3").unwrap();
    }

    #[test]
    #[should_panic]
    fn fail_non_numeric_gi() {
        let _ = parse_read_header("abc-123").unwrap();
    }

    #[test]
    #[should_panic]
    fn fail_non_numeric_taxid() {
        let _ = parse_read_header("123-abc").unwrap();
    }
}
