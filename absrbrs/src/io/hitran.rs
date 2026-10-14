//! HITRAN 160-character `.par` records.

use std::ops::Range;
use std::path::Path;

use super::{Builder, Lines, parse_error};
use crate::Result;
use crate::lines::LineList;

/// The length of a HITRAN record.
const RECORD_LEN: usize = 160;

/// Byte ranges of the fields of a record that the readers use.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct FixedLayout {
    pub(super) molecule: Range<usize>,
    pub(super) isotopologue: Range<usize>,
    pub(super) nu: Range<usize>,
    pub(super) s_ref: Range<usize>,
    pub(super) gamma_air: Range<usize>,
    pub(super) gamma_self: Range<usize>,
    pub(super) e_lower: Range<usize>,
    pub(super) n_air: Range<usize>,
    pub(super) delta_air: Range<usize>,
}

impl FixedLayout {
    /// The HITRAN 160-character layout, as in HAPI's `HITRAN_DEFAULT_HEADER`.
    pub(super) const HITRAN: FixedLayout = FixedLayout {
        molecule: 0..2,
        isotopologue: 2..3,
        nu: 3..15,
        s_ref: 15..25,
        gamma_air: 35..40,
        gamma_self: 40..45,
        e_lower: 45..55,
        n_air: 55..59,
        delta_air: 59..67,
    };

    /// The shortest record that holds every field.
    pub(super) fn min_len(&self) -> usize {
        [
            &self.molecule,
            &self.isotopologue,
            &self.nu,
            &self.s_ref,
            &self.gamma_air,
            &self.gamma_self,
            &self.e_lower,
            &self.n_air,
            &self.delta_air,
        ]
        .into_iter()
        .map(|r| r.end)
        .max()
        .unwrap_or(0)
    }
}

/// The fields of one record that the readers use, in HITRAN units.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) struct FixedRecord {
    pub(super) molecule: u16,
    pub(super) isotopologue: u16,
    pub(super) nu: f64,
    pub(super) s_ref: f64,
    pub(super) gamma_air: f64,
    pub(super) gamma_self: f64,
    pub(super) e_lower: f64,
    pub(super) n_air: f64,
    pub(super) delta_air: f64,
}

impl FixedRecord {
    /// Whether the record is kept: S_ref = 0 (as in sasktran2) and ν₀ ≤ 0 (which a line list
    /// cannot hold) are skipped.
    ///
    /// # Errors
    ///
    /// The message for [`Error::Parse`](crate::Error::Parse) if S_ref or a half-width is
    /// negative (checked for every record, kept or not).
    pub(super) fn keep(&self) -> std::result::Result<bool, String> {
        for (name, value) in [
            ("s_ref", self.s_ref),
            ("gamma_air", self.gamma_air),
            ("gamma_self", self.gamma_self),
        ] {
            if value < 0.0 {
                return Err(format!("{name} is negative ({value})"));
            }
        }
        Ok(self.s_ref != 0.0 && self.nu > 0.0)
    }
}

/// A finite number, also in the Fortran forms HAPI accepts: a `D` exponent (`1.5D-03`) and an
/// exponent with only its sign (`1.234-123` for `1.234E-123`). NaN and infinities are `None`.
pub(super) fn parse_number(text: &str) -> Option<f64> {
    let value = match text.parse::<f64>() {
        Ok(v) => v,
        Err(_) => fortran_number(text)?,
    };
    value.is_finite().then_some(value)
}

/// The Fortran forms of [`parse_number`]; only reached when Rust's syntax fails.
fn fortran_number(text: &str) -> Option<f64> {
    if text.contains(['D', 'd']) {
        return text.replace(['D', 'd'], "E").parse().ok();
    }
    if text.contains(['E', 'e']) {
        return None;
    }
    // A sign after the first character that follows a digit or a point starts the exponent.
    let at = text.get(1..)?.rfind(['+', '-'])? + 1;
    let mantissa = text.get(..at)?;
    if !mantissa.ends_with(|c: char| c.is_ascii_digit() || c == '.') {
        return None;
    }
    format!("{mantissa}E{}", text.get(at..)?).parse().ok()
}

/// The trimmed text of a field.
fn text<'a>(
    record: &'a str,
    range: &Range<usize>,
    name: &str,
) -> std::result::Result<&'a str, String> {
    record.get(range.clone()).map(str::trim).ok_or_else(|| {
        format!(
            "record is {} bytes, too short for {name} at bytes {}..{}",
            record.len(),
            range.start,
            range.end
        )
    })
}

/// A finite floating-point field.
fn float(record: &str, range: &Range<usize>, name: &str) -> std::result::Result<f64, String> {
    let t = text(record, range, name)?;
    parse_number(t).ok_or_else(|| format!("{name}: cannot parse {t:?}"))
}

/// The isotopologue number of a code: `1`–`9` are 1–9, `0` is 10, and `A`–`Z` are 11–36.
/// A code of several digits (from a HAPI header with a wider field) is the number itself.
fn isotopologue_number(code: &str) -> Option<u16> {
    match code.as_bytes() {
        [c @ b'1'..=b'9'] => Some(u16::from(c - b'0')),
        [b'0'] => Some(10),
        [c @ b'A'..=b'Z'] => Some(u16::from(c - b'A') + 11),
        [_] | [] => None,
        _ => code.parse::<u16>().ok().filter(|&n| n > 0),
    }
}

/// Parses the fields of `layout` from a record. The error is the message for
/// [`Error::Parse`](crate::Error::Parse).
pub(super) fn parse_fixed(
    record: &str,
    layout: &FixedLayout,
) -> std::result::Result<FixedRecord, String> {
    let molecule = text(record, &layout.molecule, "molecule")?;
    let molecule = molecule
        .parse::<u16>()
        .map_err(|_| format!("molecule: cannot parse {molecule:?}"))?;
    let iso = text(record, &layout.isotopologue, "isotopologue")?;
    let isotopologue =
        isotopologue_number(iso).ok_or_else(|| format!("isotopologue: cannot parse {iso:?}"))?;
    Ok(FixedRecord {
        molecule,
        isotopologue,
        nu: float(record, &layout.nu, "nu")?,
        s_ref: float(record, &layout.s_ref, "s_ref")?,
        gamma_air: float(record, &layout.gamma_air, "gamma_air")?,
        gamma_self: float(record, &layout.gamma_self, "gamma_self")?,
        e_lower: float(record, &layout.e_lower, "e_lower")?,
        n_air: float(record, &layout.n_air, "n_air")?,
        delta_air: float(record, &layout.delta_air, "delta_air")?,
    })
}

/// Reads a file of HITRAN 160-character records (`.par`).
///
/// Records with S_ref = 0 or ν₀ ≤ 0 are skipped, and so are empty lines. The optional fields
/// (n_self, δ_self, δ′_air) are absent: a `.par` record does not carry them.
///
/// # Errors
///
/// - [`Error::Io`](crate::Error::Io) if the file cannot be read;
/// - [`Error::Parse`](crate::Error::Parse), with the line number, if a record is not ASCII,
///   holds a stray carriage return, is shorter than 160 bytes, has a field that does not parse
///   (or is not finite), or has a negative S_ref or half-width.
///
/// Numbers may also take the Fortran forms `1.5D-03` and `1.234-123` (for `1.234E-123`).
pub fn read_hitran_par(path: impl AsRef<Path>) -> Result<LineList> {
    let mut lines = Lines::open(path.as_ref())?;
    let path = lines.path.clone();
    let mut builder = Builder::new([false; 3]);
    while let Some((line_number, record)) = lines.next_line()? {
        if record.len() < RECORD_LEN {
            return Err(parse_error(
                &path,
                line_number,
                format!(
                    "record is {} bytes; HITRAN records are {RECORD_LEN}",
                    record.len()
                ),
            ));
        }
        let fields = parse_fixed(record, &FixedLayout::HITRAN)
            .map_err(|message| parse_error(&path, line_number, message))?;
        let keep = fields
            .keep()
            .map_err(|message| parse_error(&path, line_number, message))?;
        if keep {
            builder.push(&fields, [f64::NAN; 3]);
        }
    }
    builder.finish()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn isotopologue_codes() {
        assert_eq!(isotopologue_number("1"), Some(1));
        assert_eq!(isotopologue_number("9"), Some(9));
        assert_eq!(isotopologue_number("0"), Some(10));
        assert_eq!(isotopologue_number("A"), Some(11));
        assert_eq!(isotopologue_number("B"), Some(12));
        assert_eq!(isotopologue_number("Z"), Some(36));
        assert_eq!(isotopologue_number("12"), Some(12));
        for bad in ["", "a", "#", "-", "00", "-1"] {
            assert_eq!(isotopologue_number(bad), None, "{bad:?}");
        }
    }

    #[test]
    fn fortran_number_forms() {
        assert_eq!(parse_number("1.5D-03"), Some(1.5e-3));
        assert_eq!(parse_number("-2.5d+2"), Some(-250.0));
        assert_eq!(parse_number("1.234-123"), Some(1.234e-123));
        assert_eq!(parse_number("1.234+05"), Some(1.234e5));
        assert_eq!(parse_number("-.5-1"), Some(-0.05));
        assert_eq!(parse_number("5.-1"), Some(0.5));
        assert_eq!(parse_number("-.004500"), Some(-0.0045));
        for bad in [
            "",
            "-",
            "+5-",
            "1.2E-3-4",
            "nan",
            "inf",
            "-infinity",
            "1D400",
            "1..2-3",
            "abc-1",
            "1-2-3",
            "D5",
            "é-1",
        ] {
            assert_eq!(parse_number(bad), None, "{bad:?}");
        }
        // Underflow to zero is finite, as with Rust's own syntax.
        assert_eq!(parse_number("1.0-400"), Some(0.0));
    }

    #[test]
    fn negative_intensities_and_widths_are_errors() {
        let r = FixedRecord {
            molecule: 1,
            isotopologue: 1,
            nu: 1.0,
            s_ref: 1e-20,
            gamma_air: 0.07,
            gamma_self: 0.09,
            e_lower: 1.0,
            n_air: 0.7,
            delta_air: 0.0,
        };
        assert_eq!(r.keep(), Ok(true));
        assert_eq!(FixedRecord { s_ref: 0.0, ..r }.keep(), Ok(false));
        assert_eq!(FixedRecord { nu: -1.0, ..r }.keep(), Ok(false));
        assert_eq!(
            FixedRecord {
                gamma_self: -0.01,
                ..r
            }
            .keep(),
            Err("gamma_self is negative (-0.01)".into())
        );
        assert!(FixedRecord { s_ref: -1e-20, ..r }.keep().is_err());
        assert!(
            FixedRecord {
                gamma_air: -0.0001,
                ..r
            }
            .keep()
            .is_err()
        );
    }

    #[test]
    fn layout_length() {
        assert_eq!(FixedLayout::HITRAN.min_len(), 67);
    }

    #[test]
    fn short_or_non_finite_fields_are_errors() {
        let record = format!("{:<67}", " 21  100.000000 1.000E-20");
        assert!(
            parse_fixed(&record[..20], &FixedLayout::HITRAN)
                .unwrap_err()
                .contains("too short for s_ref")
        );
        let mut nan = format!("{:<160}", " 21  100.000000       nan");
        assert_eq!(
            parse_fixed(&nan, &FixedLayout::HITRAN).unwrap_err(),
            "s_ref: cannot parse \"nan\""
        );
        nan.replace_range(15..25, " 1.000E-20");
        assert_eq!(
            parse_fixed(&nan, &FixedLayout::HITRAN).unwrap_err(),
            "gamma_air: cannot parse \"\""
        );
    }
}
