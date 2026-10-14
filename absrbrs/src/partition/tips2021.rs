//! The built-in TIPS-2021 tables, embedded from `data/tips2021.bin` (written by
//! `tools/tips2021.py`; format below).
//!
//! All little-endian: the 16-byte magic `b"ABSRBRS-TIPS2021"`; a `u32` count; `count` entries of
//! 16 bytes (`u16` molecule, `u16` isotopologue, `u32` n_points, `f32` t_first, `f32` t_step);
//! then the `f32` Q values of each entry in turn.

use std::collections::BTreeMap;
use std::fmt;
use std::sync::OnceLock;

use super::hermite::Table;
use super::{PartitionFunction, PartitionValue, check_range};
use crate::{Error, Result};

static BLOB: &[u8] = include_bytes!("../../data/tips2021.bin");
const MAGIC: &[u8; 16] = b"ABSRBRS-TIPS2021";

/// One isotopologue's table: where its values are in `BLOB`, and the grid.
#[derive(Debug)]
struct Entry {
    /// Byte offset of the first Q value in `BLOB`.
    offset: usize,
    n: usize,
    t_first: f64,
    t_step: f64,
    /// Decoded on first use.
    table: OnceLock<Table>,
}

impl Entry {
    fn t_last(&self) -> f64 {
        self.t_first + (self.n - 1) as f64 * self.t_step
    }

    fn table(&self) -> &Table {
        self.table.get_or_init(|| {
            let bytes = &BLOB[self.offset..self.offset + 4 * self.n];
            let t = (0..self.n)
                .map(|k| self.t_first + k as f64 * self.t_step)
                .collect();
            let ln_q = bytes
                .chunks_exact(4)
                .map(|b| f64::from(f32::from_le_bytes(le4(b))).ln())
                .collect();
            Table::new(t, ln_q)
        })
    }
}

/// Total internal partition sums from TIPS-2021 (Gamache et al. 2021, doi:10.1016/j.jqsrt.2021.107713;
/// MIT licence, Copyright (C) 2021 Robert Gamache), for 181 HITRAN isotopologues at 1 K steps
/// from 1 K to each table's maximum (1000–9000 K). Q and d ln Q/dT come from cubic Hermite
/// interpolation of ln Q in T, so they are consistent with each other.
///
/// Three tables start above 1 K. The H2S tables (molecule 31, isotopologues 2 and 3) start at
/// 30 K: below that their ratios to isotopologue 1 depart from the ratio of the degeneracies, by
/// up to a factor of 2. CH3 (57, 1) starts at 2 K, because its Q at 1 K is 0, and its
/// low-temperature values carry few digits (Q(2 K) = 0.000416, about 1e-3 relative rounding).
/// `TIPS_2021.py` gives 3010 K as the maximum for CF4 (42, 1), but its table runs, smoothly, to
/// 5000 K, and all of it is used.
///
/// Accuracy:
/// - Q is stored as f32, which rounds by at most 6e-8 relative, comparable to the 6–9
///   significant figures of the source.
/// - Between the nodes, from 50 K up, the interpolant matches higher-order interpolation to about
///   1e-8 in ln Q; d ln Q/dT carries about 1e-5 relative noise from the data.
/// - Below about 20 K, where ln Q curves sharply, errors between the nodes are 1e-2 to 2e-1.
/// - Kinks in the source tables (SO3 (47, 1) at 650 K, HCN (23, 1) and (23, 2) at 1000 K)
///   distort d ln Q/dT within ±1 K of them.
/// - At non-integer temperatures the values differ from those of TIPS's own (and HAPI's)
///   interpolation, which is linear in Q, by up to about 4e-6 relative.
pub struct Tips2021 {
    index: BTreeMap<(u16, u16), Entry>,
}

impl fmt::Debug for Tips2021 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Tips2021")
            .field("isotopologues", &self.index.len())
            .finish_non_exhaustive()
    }
}

/// The built-in TIPS-2021 partition functions, parsed on first use.
pub fn tips2021() -> &'static Tips2021 {
    static TIPS: OnceLock<Tips2021> = OnceLock::new();
    TIPS.get_or_init(Tips2021::parse)
}

fn le4(b: &[u8]) -> [u8; 4] {
    b.try_into().expect("embedded tips2021.bin is malformed")
}

impl Tips2021 {
    /// Parses the header of the embedded file. The file is compile-time data whose layout
    /// `tests/test_data_files.py` checks against its manifest, so a malformed file is a build
    /// defect, not a runtime condition: the `expect`s here cannot fire on a correct build.
    fn parse() -> Tips2021 {
        const MALFORMED: &str = "embedded tips2021.bin is malformed";
        assert!(BLOB.len() >= 20 && &BLOB[..16] == MAGIC, "{MALFORMED}");
        let count = u32::from_le_bytes(le4(&BLOB[16..20])) as usize;
        let mut offset = 20 + 16 * count;
        let mut index = BTreeMap::new();
        for k in 0..count {
            let e = BLOB.get(20 + 16 * k..36 + 16 * k).expect(MALFORMED);
            let molecule = u16::from_le_bytes([e[0], e[1]]);
            let isotopologue = u16::from_le_bytes([e[2], e[3]]);
            let n = u32::from_le_bytes(le4(&e[4..8])) as usize;
            let t_first = f64::from(f32::from_le_bytes(le4(&e[8..12])));
            let t_step = f64::from(f32::from_le_bytes(le4(&e[12..16])));
            assert!(n >= 2 && t_step > 0.0, "{MALFORMED}");
            index.insert(
                (molecule, isotopologue),
                Entry {
                    offset,
                    n,
                    t_first,
                    t_step,
                    table: OnceLock::new(),
                },
            );
            offset += 4 * n;
        }
        assert!(offset == BLOB.len() && index.len() == count, "{MALFORMED}");
        Tips2021 { index }
    }

    /// The first and last temperatures (K) of an isotopologue's table, or `None` if there is
    /// none.
    pub fn temperature_range(&self, molecule: u16, isotopologue: u16) -> Option<(f64, f64)> {
        self.index
            .get(&(molecule, isotopologue))
            .map(|e| (e.t_first, e.t_last()))
    }

    /// The (molecule, isotopologue) pairs with a table, sorted.
    pub fn isotopologues(&self) -> Vec<(u16, u16)> {
        self.index.keys().copied().collect()
    }
}

impl PartitionFunction for Tips2021 {
    fn evaluate(
        &self,
        molecule: u16,
        isotopologue: u16,
        temperature: f64,
    ) -> Result<PartitionValue> {
        let entry =
            self.index
                .get(&(molecule, isotopologue))
                .ok_or(Error::MissingPartitionFunction {
                    molecule,
                    isotopologue,
                })?;
        check_range(
            (entry.t_first, entry.t_last()),
            molecule,
            isotopologue,
            temperature,
        )?;
        let (ln_q, dlnq_dt) = entry.table().evaluate(temperature);
        Ok(PartitionValue {
            q: ln_q.exp(),
            dlnq_dt,
        })
    }
}
