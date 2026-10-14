//! Error type shared by every absrbrs computation (spec §7.4).
//!
//! [`Error`] is `Clone` and `PartialEq` so that results can be compared in tests and stored
//! freely. It therefore holds only owned, comparable data: an I/O failure keeps its
//! [`std::io::ErrorKind`] and message rather than the `std::io::Error` itself, and no variant
//! carries a `source()` chain.

/// Errors returned by absrbrs.
///
/// A new variant also needs an arm in `absrbrs-py/src/errors.rs` (`to_pyerr`), a class in
/// `python/absrbrs/errors.py` if it should not surface as a plain `AbsrbrsError`, and a case in
/// the `_testing.raise_core_error` hook with a test in `tests/test_errors.py`.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
#[non_exhaustive]
pub enum Error {
    /// An input is outside its valid range: a wavenumber grid that is not strictly increasing,
    /// a non-positive temperature, a NaN, a buffer of the wrong size, an output too large, ...
    #[error("invalid input: {0}")]
    InvalidInput(String),
    /// No partition function is available for an isotopologue.
    #[error("no partition function for molecule {molecule}, isotopologue {isotopologue}")]
    MissingPartitionFunction {
        /// HITRAN molecule number.
        molecule: u16,
        /// HITRAN isotopologue number (codes 0, A and B are 10, 11 and 12).
        isotopologue: u16,
    },
    /// A temperature is outside the range of a partition-function table.
    #[error(
        "temperature {temperature:?} K is outside the partition-function table of molecule \
         {molecule}, isotopologue {isotopologue}"
    )]
    OutsideTable {
        /// HITRAN molecule number.
        molecule: u16,
        /// HITRAN isotopologue number (codes 0, A and B are 10, 11 and 12).
        isotopologue: u16,
        /// The temperature requested, in K.
        temperature: f64,
    },
    /// A record of a line file is malformed.
    #[error("{path}, line {line_number}: {message}")]
    Parse {
        /// The file, as given, or a description such as `<string>` for input from memory.
        path: String,
        /// One-based line number of the record in its file.
        line_number: usize,
        /// What is wrong with the record.
        message: String,
    },
    /// A file could not be read.
    #[error("cannot read {path}: {message}")]
    Io {
        /// The path, as given.
        path: String,
        /// The kind of failure, which tells a missing file from other failures.
        kind: std::io::ErrorKind,
        /// The operating system's description of the failure.
        message: String,
    },
}

/// Result alias using [`Error`].
pub type Result<T> = core::result::Result<T, Error>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn messages_are_informative() {
        assert_eq!(
            Error::InvalidInput("grid is not increasing".into()).to_string(),
            "invalid input: grid is not increasing"
        );
        assert_eq!(
            Error::MissingPartitionFunction {
                molecule: 2,
                isotopologue: 11
            }
            .to_string(),
            "no partition function for molecule 2, isotopologue 11"
        );
        assert_eq!(
            Error::OutsideTable {
                molecule: 2,
                isotopologue: 1,
                temperature: 6000.5
            }
            .to_string(),
            "temperature 6000.5 K is outside the partition-function table of molecule 2, \
             isotopologue 1"
        );
        assert_eq!(
            Error::Parse {
                path: "CO2.par".into(),
                line_number: 17,
                message: "bad width".into()
            }
            .to_string(),
            "CO2.par, line 17: bad width"
        );
        assert_eq!(
            Error::Io {
                path: "CO2.par".into(),
                kind: std::io::ErrorKind::NotFound,
                message: "No such file".into()
            }
            .to_string(),
            "cannot read CO2.par: No such file"
        );
    }

    #[test]
    fn extreme_temperatures_print_compactly() {
        let message = Error::OutsideTable {
            molecule: 1,
            isotopologue: 1,
            temperature: 1e300,
        }
        .to_string();
        assert!(message.starts_with("temperature 1e300 K "), "{message}");
    }

    #[test]
    fn error_is_send_sync_and_static() {
        fn check<T: Send + Sync + 'static>() {}
        check::<Error>();
    }
}
