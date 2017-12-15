//! Result and Error types for all vedro code.

use bincode::rustc_serialize::{DecodingError, EncodingError};
use std::fmt;
use std::io;
use std::str;

#[allow(missing_docs)]
pub type VedroResult<T> = Result<T, VedroError>;

#[allow(missing_docs)]
#[derive(Debug)]
pub enum VedroError {
    Io(io::Error),
    Deserialize(DecodingError),
    InvalidHeader(String),
    InvalidInteger(String),
    MissingFile(String),
    MissingHeader,
    Serialize(EncodingError),
    Utf8(str::Utf8Error),
}

impl fmt::Display for VedroError {
    fn fmt(&self, f: &mut fmt::Formatter) -> Result<(), fmt::Error> {

        match self {
            &VedroError::Io(ref e) => write!(f, "I/O problem: {}", e),
            &VedroError::Deserialize(ref e) => {
                write!(f, "Unable to deserialize item from disk: {}", e)
            },
            &VedroError::InvalidHeader(ref h) => {
                write!(f, "Incorrectly formatted FASTA header: {}", h)
            },
            &VedroError::InvalidInteger(ref s) => write!(f, "Unable to parse \"{}\" as integer", s),
            &VedroError::MissingFile(ref p) => write!(f, "Unable to find file {}", p),
            &VedroError::MissingHeader => write!(f, "Empty header found in FASTA file"),
            &VedroError::Serialize(ref e) => write!(f, "Unable to serialize item to disk: {}", e),
            &VedroError::Utf8(ref e) => write!(f, "Found invalid UTF8 input ({})", e),
        }
    }
}

impl From<io::Error> for VedroError {
    fn from(e: io::Error) -> Self {
        VedroError::Io(e)
    }
}

impl From<DecodingError> for VedroError {
    fn from(e: DecodingError) -> Self {
        VedroError::Deserialize(e)
    }
}

impl From<EncodingError> for VedroError {
    fn from(e: EncodingError) -> Self {
        VedroError::Serialize(e)
    }
}

impl From<str::Utf8Error> for VedroError {
    fn from(e: str::Utf8Error) -> Self {
        VedroError::Utf8(e)
    }
}
