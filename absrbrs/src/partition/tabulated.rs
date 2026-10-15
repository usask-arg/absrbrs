//! User-supplied partition-function tables.

use std::collections::BTreeMap;

use super::hermite::Table;
use super::{PartitionFunction, PartitionValue, check_range};
use crate::{Error, Result};

/// Partition functions from user-supplied tables of Q(T), by cubic Hermite interpolation of
/// ln Q against ln T (exact for Q ∝ Tⁿ, which suits coarse tables); d ln Q/dT is
/// (d ln Q/d ln T)/T. The temperature grid need not be uniform. ([`Tips2021`], with its 1 K
/// grid, interpolates in T instead.)
///
/// [`Tips2021`]: super::Tips2021
#[derive(Clone, Debug, Default)]
pub struct TabulatedPartition {
    tables: BTreeMap<(u16, u16), UserTable>,
}

/// A table in ln T, with its temperature range kept exactly as given.
#[derive(Clone, Debug)]
struct UserTable {
    range: (f64, f64),
    table: Table,
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
        let ln_t: Vec<f64> = temperature.iter().map(|t| t.ln()).collect();
        if !ln_t.windows(2).all(|w| w[0] < w[1]) {
            return invalid(
                "temperatures are too close together (ln T must be strictly increasing)",
            );
        }
        let ln_q = q.iter().map(|q| q.ln()).collect();
        let range = (temperature[0], temperature[temperature.len() - 1]);
        let table = UserTable {
            range,
            table: Table::new(ln_t, ln_q),
        };
        self.tables.insert((molecule, isotopologue), table);
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
        check_range(table.range, molecule, isotopologue, temperature)?;
        let (ln_q, dlnq_dlnt) = table.table.evaluate(temperature.ln());
        Ok(PartitionValue {
            q: ln_q.exp(),
            dlnq_dt: dlnq_dlnt / temperature,
        })
    }
}
