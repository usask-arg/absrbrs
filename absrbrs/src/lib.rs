//! Molecular absorption cross sections with first-class linearization.
//!
//! The design is described in `docs/specs/2026-10-14-absrbrs-design.md` in the repository.
#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod faddeeva;
#[cfg(test)]
mod test_data;

pub use error::{Error, Result};
