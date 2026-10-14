//! Error type shared by every absrbrs computation (spec §7.4).

/// Errors returned by absrbrs.
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
        "temperature {temperature} K is outside the partition-function table of molecule \
         {molecule}, isotopologue {isotopologue}"
    )]
    OutsideTable {
        /// HITRAN molecule number.
        molecule: u16,
        /// HITRAN isotopologue number.
        isotopologue: u16,
        /// The temperature requested, in K.
        temperature: f64,
    },
    /// A record of a line file is malformed.
    #[error("line {line_number}: {message}")]
    Parse {
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
                line_number: 17,
                message: "bad width".into()
            }
            .to_string(),
            "line 17: bad width"
        );
        assert_eq!(
            Error::Io {
                path: "CO2.par".into(),
                message: "No such file".into()
            }
            .to_string(),
            "cannot read CO2.par: No such file"
        );
    }
}
