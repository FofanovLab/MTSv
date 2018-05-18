//! Result and Error types for all mtsv code.

use bincode::rustc_serialize::{DecodingError, EncodingError};
use std::fmt;
use std::io;
use std::str;

#[allow(missing_docs)]
pub type MtsvResult<T> = Result<T, MtsvError>;

#[allow(missing_docs)]
#[derive(Debug)]
pub enum MtsvError {
    Io(io::Error),
    Deserialize(DecodingError),
    InvalidHeader(String),
    InvalidInteger(String),
    MissingFile(String),
    MissingHeader,
    Serialize(EncodingError),
    Utf8(str::Utf8Error),
}

impl fmt::Display for MtsvError {
    fn fmt(&self, f: &mut fmt::Formatter) -> Result<(), fmt::Error> {

        match self {
            &MtsvError::Io(ref e) => write!(f, "I/O problem: {}", e),
            &MtsvError::Deserialize(ref e) => {
                write!(f, "Unable to deserialize item from disk: {}", e)
            },
            &MtsvError::InvalidHeader(ref h) => {
                write!(f, "Incorrectly formatted FASTA header: {}", h)
            },
            &MtsvError::InvalidInteger(ref s) => write!(f, "Unable to parse \"{}\" as integer", s),
            &MtsvError::MissingFile(ref p) => write!(f, "Unable to find file {}", p),
            &MtsvError::MissingHeader => write!(f, "Empty header found in FASTA file"),
            &MtsvError::Serialize(ref e) => write!(f, "Unable to serialize item to disk: {}", e),
            &MtsvError::Utf8(ref e) => write!(f, "Found invalid UTF8 input ({})", e),
        }
    }
}

impl From<io::Error> for MtsvError {
    fn from(e: io::Error) -> Self {
        MtsvError::Io(e)
    }
}

impl From<DecodingError> for MtsvError {
    fn from(e: DecodingError) -> Self {
        MtsvError::Deserialize(e)
    }
}

impl From<EncodingError> for MtsvError {
    fn from(e: EncodingError) -> Self {
        MtsvError::Serialize(e)
    }
}

impl From<str::Utf8Error> for MtsvError {
    fn from(e: str::Utf8Error) -> Self {
        MtsvError::Utf8(e)
    }
}
