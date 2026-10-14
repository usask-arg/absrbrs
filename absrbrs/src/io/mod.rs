//! Line-file readers: HITRAN `.par` records and HAPI `.data`/`.header` caches (spec §5.6).
//!
//! Both readers skip records with S_ref = 0 (as sasktran2 does) and records with ν₀ ≤ 0, which
//! a [`LineList`] cannot hold.

mod hapi;
mod hitran;

pub use hapi::read_hapi;
pub use hitran::read_hitran_par;

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use crate::lines::{LineColumns, LineList};
use crate::{Error, Result};
use hitran::FixedRecord;

/// The [`Error::Io`] for a failure on `path`.
fn io_error(path: &str, e: &std::io::Error) -> Error {
    Error::Io {
        path: path.to_owned(),
        kind: e.kind(),
        message: e.to_string(),
    }
}

/// Reads the non-blank lines of a text file one at a time into a reused buffer.
struct Lines {
    reader: BufReader<File>,
    buffer: Vec<u8>,
    line_number: usize,
    path: String,
}

impl Lines {
    fn open(path: &Path) -> Result<Lines> {
        let shown = path.display().to_string();
        let file = File::open(path).map_err(|e| io_error(&shown, &e))?;
        Ok(Lines {
            reader: BufReader::with_capacity(1 << 16, file),
            buffer: Vec::with_capacity(512),
            line_number: 0,
            path: shown,
        })
    }

    /// The next non-blank line without its line ending (`\n` or `\r\n`), and its one-based
    /// line number.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] if reading fails, and [`Error::Parse`] if the line is not ASCII or holds
    /// a carriage return that is not part of its line ending.
    fn next_line(&mut self) -> Result<Option<(usize, &str)>> {
        let len = loop {
            self.buffer.clear();
            let n = self
                .reader
                .read_until(b'\n', &mut self.buffer)
                .map_err(|e| io_error(&self.path, &e))?;
            if n == 0 {
                return Ok(None);
            }
            self.line_number += 1;
            let mut line = self.buffer.as_slice();
            if let Some(l) = line.strip_suffix(b"\n") {
                line = l.strip_suffix(b"\r").unwrap_or(l);
            }
            // Lines that are blank (or a lone `\r` at the end of the file) are skipped.
            if !line.trim_ascii_end().is_empty() {
                break line.len();
            }
        };
        let line = &self.buffer[..len];
        if line.contains(&b'\r') {
            // Old Mac (CR-only) line endings would otherwise hide records inside one line.
            return Err(parse_error(
                &self.path,
                self.line_number,
                "contains a carriage return inside the line".into(),
            ));
        }
        if !line.is_ascii() {
            return Err(parse_error(
                &self.path,
                self.line_number,
                "contains non-ASCII bytes".into(),
            ));
        }
        // ASCII is valid UTF-8, so this cannot fail.
        let line = std::str::from_utf8(line)
            .map_err(|e| parse_error(&self.path, self.line_number, e.to_string()))?;
        Ok(Some((self.line_number, line)))
    }
}

/// The [`Error::Parse`] for line `line_number` of `path`.
fn parse_error(path: &str, line_number: usize, message: String) -> Error {
    Error::Parse {
        path: path.to_owned(),
        line_number,
        message,
    }
}

/// Collects kept records into [`LineColumns`].
struct Builder {
    columns: LineColumns,
}

impl Builder {
    /// `optional` says which of n_self, δ_self and δ′_air have columns.
    fn new(optional: [bool; 3]) -> Builder {
        let column = |present: bool| present.then(Vec::new);
        Builder {
            columns: LineColumns {
                n_self: column(optional[0]),
                delta_self: column(optional[1]),
                deltap_air: column(optional[2]),
                ..LineColumns::default()
            },
        }
    }

    /// Appends a record and its n_self, δ_self and δ′_air (NaN when absent; ignored for a
    /// field without a column).
    fn push(&mut self, r: &FixedRecord, optional: [f64; 3]) {
        let c = &mut self.columns;
        c.molecule.push(r.molecule);
        c.isotopologue.push(r.isotopologue);
        c.nu.push(r.nu);
        c.s_ref.push(r.s_ref);
        c.e_lower.push(r.e_lower);
        c.gamma_air.push(r.gamma_air);
        c.gamma_self.push(r.gamma_self);
        c.n_air.push(r.n_air);
        c.delta_air.push(r.delta_air);
        for (column, value) in [&mut c.n_self, &mut c.delta_self, &mut c.deltap_air]
            .into_iter()
            .zip(optional)
        {
            if let Some(column) = column {
                column.push(value);
            }
        }
    }

    fn finish(self) -> Result<LineList> {
        LineList::from_columns(self.columns)
    }
}
