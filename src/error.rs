//! Result and Error types for all mtsv code.

use bincode::rustc_serialize::{DecodingError, EncodingError};
use std::fmt;
use std::io;
use std::str;

#[allow(missing_docs)]
pub type mtsvResult<T> = Result<T, mtsvError>;

#[allow(missing_docs)]
#[derive(Debug)]
pub enum mtsvError {
    Io(io::Error),
    Deserialize(DecodingError),
    InvalidHeader(String),
    InvalidInteger(String),
    MissingFile(String),
    MissingHeader,
    Serialize(EncodingError),
    Utf8(str::Utf8Error),
}

impl fmt::Display for mtsvError {
    fn fmt(&self, f: &mut fmt::Formatter) -> Result<(), fmt::Error> {

        match self {
            &mtsvError::Io(ref e) => write!(f, "I/O problem: {}", e),
            &mtsvError::Deserialize(ref e) => {
                write!(f, "Unable to deserialize item from disk: {}", e)
            },
            &mtsvError::InvalidHeader(ref h) => {
                write!(f, "Incorrectly formatted FASTA header: {}", h)
            },
            &mtsvError::InvalidInteger(ref s) => write!(f, "Unable to parse \"{}\" as integer", s),
            &mtsvError::MissingFile(ref p) => write!(f, "Unable to find file {}", p),
            &mtsvError::MissingHeader => write!(f, "Empty header found in FASTA file"),
            &mtsvError::Serialize(ref e) => write!(f, "Unable to serialize item to disk: {}", e),
            &mtsvError::Utf8(ref e) => write!(f, "Found invalid UTF8 input ({})", e),
        }
    }
}

impl From<io::Error> for mtsvError {
    fn from(e: io::Error) -> Self {
        mtsvError::Io(e)
    }
}

impl From<DecodingError> for mtsvError {
    fn from(e: DecodingError) -> Self {
        mtsvError::Deserialize(e)
    }
}

impl From<EncodingError> for mtsvError {
    fn from(e: EncodingError) -> Self {
        mtsvError::Serialize(e)
    }
}

impl From<str::Utf8Error> for mtsvError {
    fn from(e: str::Utf8Error) -> Self {
        mtsvError::Utf8(e)
    }
}
