//! User-supplied partition-function tables.

use std::collections::BTreeMap;

use super::hermite::Table;
use super::{PartitionFunction, PartitionValue, evaluate_table};
use crate::{Error, Result};

/// Partition functions from user-supplied tables of Q(T), interpolated like [`Tips2021`]
/// (cubic Hermite in ln Q). The temperature grid need not be uniform.
///
/// [`Tips2021`]: super::Tips2021
#[derive(Clone, Debug, Default)]
pub struct TabulatedPartition {
    tables: BTreeMap<(u16, u16), Table>,
}

impl TabulatedPartition {
    /// An empty set of tables.
    pub fn new() -> TabulatedPartition {
        TabulatedPartition::default()
    }

    /// Adds the table of Q at `temperature` (K) for an isotopologue, replacing any previous
    /// table for it.
    ///
    /// Errors with `InvalidInput` unless `temperature` and `q` have the same length of at least
    /// two, `temperature` is finite, positive and strictly increasing, and `q` is finite and
    /// positive.
    pub fn insert(
        &mut self,
        molecule: u16,
        isotopologue: u16,
        temperature: &[f64],
        q: &[f64],
    ) -> Result<()> {
        let invalid = |what: &str| {
            Err(Error::InvalidInput(format!(
                "partition table for molecule {molecule}, isotopologue {isotopologue}: {what}"
            )))
        };
        if temperature.len() != q.len() {
            return invalid(&format!(
                "{} temperatures but {} values of Q",
                temperature.len(),
                q.len()
            ));
        }
        if temperature.len() < 2 {
            return invalid("at least two points are needed");
        }
        if !temperature.iter().all(|t| t.is_finite() && *t > 0.0) {
            return invalid("temperatures must be finite and positive");
        }
        if !temperature.windows(2).all(|w| w[0] < w[1]) {
            return invalid("temperatures must be strictly increasing");
        }
        if !q.iter().all(|q| q.is_finite() && *q > 0.0) {
            return invalid("Q must be finite and positive");
        }
        let ln_q = q.iter().map(|q| q.ln()).collect();
        self.tables.insert(
            (molecule, isotopologue),
            Table::new(temperature.to_vec(), ln_q),
        );
        Ok(())
    }
}

impl PartitionFunction for TabulatedPartition {
    fn evaluate(
        &self,
        molecule: u16,
        isotopologue: u16,
        temperature: f64,
    ) -> Result<PartitionValue> {
        let table =
            self.tables
                .get(&(molecule, isotopologue))
                .ok_or(Error::MissingPartitionFunction {
                    molecule,
                    isotopologue,
                })?;
        evaluate_table(table, molecule, isotopologue, temperature)
    }
}
