//! HAPI caches: a JSON `.header` describing a `.data` file of fixed-width records followed by
//! separated extra columns.

use std::collections::HashMap;
use std::ops::Range;
use std::path::Path;

use serde_json::{Map, Value};

use super::hitran::{FixedLayout, parse_fixed, parse_number};
use super::{Builder, Lines, io_error, parse_error};
use crate::Result;
use crate::lines::LineList;

/// The optional fields, in [`Builder`] order.
const OPTIONAL: [&str; 3] = ["n_self", "delta_self", "deltap_air"];

/// Where a record holds an optional field.
#[derive(Clone, Debug, PartialEq, Eq)]
enum Source {
    /// A fixed-width column (in `order`, as `hapi.cache2storage` writes them).
    Fixed(Range<usize>),
    /// The index of an extra column (in `extra`, as `hapi.fetch` writes them).
    Extra(usize),
}

/// What the reader needs from a `.header`.
#[derive(Debug)]
struct Header {
    /// Byte ranges of the fixed fields that are read.
    layout: FixedLayout,
    /// Where the fixed-width part of a record ends (the extras follow the next separator).
    fixed_end: usize,
    separator: String,
    /// Where each of [`OPTIONAL`] is, if anywhere.
    optional: [Option<Source>; 3],
}

/// The width of a printf-style format such as `%12.6f`.
fn format_width(format: &str) -> Option<usize> {
    let rest = format.strip_prefix('%')?;
    let rest = rest.trim_start_matches(['-', '+', ' ', '#']);
    let digits = rest.len() - rest.trim_start_matches(|c: char| c.is_ascii_digit()).len();
    rest[..digits].parse::<usize>().ok().filter(|&w| w > 0)
}

impl Header {
    /// Interprets a header; the error is the message for [`Error::Parse`](crate::Error::Parse).
    fn parse(text: &str) -> std::result::Result<Header, String> {
        let value: Value =
            serde_json::from_str(text).map_err(|e| format!("invalid JSON header: {e}"))?;
        let strings = |key: &str, required: bool| -> std::result::Result<Vec<&str>, String> {
            match value.get(key) {
                None if !required => Ok(Vec::new()),
                Some(Value::Array(items)) => items
                    .iter()
                    .map(|v| v.as_str().ok_or(format!("`{key}` holds a non-string")))
                    .collect(),
                _ => Err(format!("the header has no `{key}` list")),
            }
        };
        let order = strings("order", true)?;
        let extra = strings("extra", false)?;
        let empty = Map::new();
        let format = value
            .get("format")
            .and_then(Value::as_object)
            .ok_or("the header has no `format` object")?;
        let position = value
            .get("position")
            .map(|p| p.as_object().unwrap_or(&empty));

        // The byte range of each field that has a width and a start.
        let mut spans: HashMap<&str, Range<usize>> = HashMap::new();
        let mut running = Some(0usize);
        for &name in &order {
            let width = format
                .get(name)
                .and_then(Value::as_str)
                .and_then(format_width);
            let start = match position {
                Some(p) => p
                    .get(name)
                    .and_then(Value::as_u64)
                    .and_then(|s| usize::try_from(s).ok()),
                None => running,
            };
            let end = start.zip(width).and_then(|(s, w)| s.checked_add(w));
            running = end;
            if let (Some(start), Some(end)) = (start, end) {
                spans.insert(name, start..end);
            }
        }
        let span = |name: &str| {
            spans
                .get(name)
                .cloned()
                .ok_or_else(|| format!("the header gives no fixed-width field `{name}`"))
        };
        let layout = FixedLayout {
            molecule: span("molec_id")?,
            isotopologue: span("local_iso_id")?,
            nu: span("nu")?,
            s_ref: span("sw")?,
            gamma_air: span("gamma_air")?,
            gamma_self: span("gamma_self")?,
            e_lower: span("elower")?,
            n_air: span("n_air")?,
            delta_air: span("delta_air")?,
        };
        let fixed_end = spans.values().map(|r| r.end).max().unwrap_or(0);
        let separator = match value.get("extra_separator") {
            None => ",".to_owned(),
            Some(Value::String(s)) if !s.is_empty() => s.clone(),
            Some(_) => return Err("`extra_separator` is not a non-empty string".into()),
        };
        // A field in both `order` and `extra` is read from its fixed-width column.
        let optional = OPTIONAL.map(|name| {
            spans
                .get(name)
                .cloned()
                .map(Source::Fixed)
                .or_else(|| extra.iter().position(|&e| e == name).map(Source::Extra))
        });
        Ok(Header {
            layout,
            fixed_end,
            separator,
            optional,
        })
    }

    /// The n_self, δ_self and δ′_air of a record: NaN when absent from the header, missing from
    /// the record, or not a finite number (HAPI writes `#` for missing).
    fn optional_values(&self, record: &str) -> [f64; 3] {
        let number = |token: &str| parse_number(token.trim()).unwrap_or(f64::NAN);
        let mut values = [f64::NAN; 3];
        let mut last_extra = None;
        for (value, source) in values.iter_mut().zip(&self.optional) {
            match source {
                Some(Source::Fixed(range)) => {
                    *value = record.get(range.clone()).map_or(f64::NAN, number);
                }
                Some(Source::Extra(k)) => last_extra = last_extra.max(Some(*k)),
                None => {}
            }
        }
        let Some(last) = last_extra else {
            return values;
        };
        // The extras start after the first separator that follows the fixed-width part.
        let tail = record.get(self.fixed_end.min(record.len())..).unwrap_or("");
        let Some(at) = tail.find(self.separator.as_str()) else {
            return values;
        };
        let extras = &tail[at + self.separator.len()..];
        for (k, token) in extras
            .split(self.separator.as_str())
            .take(last + 1)
            .enumerate()
        {
            for (value, source) in values.iter_mut().zip(&self.optional) {
                if *source == Some(Source::Extra(k)) {
                    *value = number(token);
                }
            }
        }
        values
    }
}

/// Reads a HAPI cache: `<directory>/<table>.header` (JSON) and `<directory>/<table>.data`.
///
/// The fixed-width fields are located with the header's `order`, `format` and `position`
/// (or, without `position`, by the running sum of the widths in `order`). n_self, δ_self and
/// δ′_air are read from fixed-width columns when `order` lists them (as `hapi.cache2storage`
/// writes them), and otherwise from the extra columns when `extra` lists them (as `hapi.fetch`
/// writes them). A missing value (`#`, or anything that is not a finite number) is NaN, meaning
/// absent for that line. A field in neither list has no column.
///
/// Records with S_ref = 0 or ν₀ ≤ 0 are skipped, and so are blank lines. Numbers may also take
/// the Fortran forms `1.5D-03` and `1.234-123` (for `1.234E-123`).
///
/// # Errors
///
/// - [`Error::Io`](crate::Error::Io) if either file cannot be read;
/// - [`Error::Parse`](crate::Error::Parse) for a header that is not UTF-8 JSON or lacks a
///   field (at line 1 of the header, or the line of an invalid UTF-8 byte), and, with the line
///   number, for a data record that is not ASCII, holds a stray carriage return, is too short
///   for its fixed fields, has a required fixed field that does not parse, or has a negative
///   S_ref or half-width.
pub fn read_hapi(directory: impl AsRef<Path>, table: &str) -> Result<LineList> {
    let directory = directory.as_ref();
    let header_path = directory.join(format!("{table}.header"));
    let shown = header_path.display().to_string();
    let bytes = std::fs::read(&header_path).map_err(|e| io_error(&shown, &e))?;
    let text = String::from_utf8(bytes).map_err(|e| {
        let valid = &e.as_bytes()[..e.utf8_error().valid_up_to()];
        let line_number = 1 + valid.iter().filter(|&&b| b == b'\n').count();
        parse_error(
            &shown,
            line_number,
            format!("header is not UTF-8: {}", e.utf8_error()),
        )
    })?;
    let header = Header::parse(&text).map_err(|message| parse_error(&shown, 1, message))?;

    let mut lines = Lines::open(&directory.join(format!("{table}.data")))?;
    let path = lines.path.clone();
    let min_len = header.layout.min_len();
    let mut builder = Builder::new(header.optional.each_ref().map(Option::is_some));
    while let Some((line_number, record)) = lines.next_line()? {
        if record.len() < min_len {
            return Err(parse_error(
                &path,
                line_number,
                format!(
                    "record is {} bytes; the header's fixed fields need {min_len}",
                    record.len()
                ),
            ));
        }
        let fields = parse_fixed(record, &header.layout)
            .map_err(|message| parse_error(&path, line_number, message))?;
        let keep = fields
            .keep()
            .map_err(|message| parse_error(&path, line_number, message))?;
        if keep {
            builder.push(&fields, header.optional_values(record));
        }
    }
    builder.finish()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn widths_of_formats() {
        assert_eq!(format_width("%12.6f"), Some(12));
        assert_eq!(format_width("%1d"), Some(1));
        assert_eq!(format_width("%-15s"), Some(15));
        assert_eq!(format_width("%s"), None);
        assert_eq!(format_width("12.6f"), None);
        assert_eq!(format_width("%0d"), None);
    }

    const ORDER: &str = r#""order": ["molec_id", "local_iso_id", "nu", "sw", "a", "gamma_air",
        "gamma_self", "elower", "n_air", "delta_air", "quanta"],
        "format": {"molec_id": "%2d", "local_iso_id": "%1d", "nu": "%12.6f", "sw": "%10.3E",
        "a": "%10.3E", "gamma_air": "%5.4f", "gamma_self": "%5.3f", "elower": "%10.4f",
        "n_air": "%4.2f", "delta_air": "%8.6f", "quanta": "%93s"}"#;

    #[test]
    fn positions_default_to_running_widths() {
        let h = Header::parse(&format!("{{{ORDER}}}")).unwrap();
        assert_eq!(h.layout, FixedLayout::HITRAN);
        assert_eq!(h.fixed_end, 160);
        assert_eq!(h.separator, ",");
        assert_eq!(h.optional, [None, None, None]);
    }

    #[test]
    fn header_errors_name_the_problem() {
        let missing = ORDER.replace("\"elower\", ", "");
        let err = Header::parse(&format!("{{{missing}}}")).unwrap_err();
        assert!(err.contains("`elower`"), "{err}");
        assert!(Header::parse("{}").unwrap_err().contains("`order`"));
        assert!(Header::parse("[1, 2").unwrap_err().contains("invalid JSON"));
        let err = Header::parse(&format!("{{{ORDER}, \"extra_separator\": \"\"}}")).unwrap_err();
        assert!(err.contains("extra_separator"), "{err}");
        // A huge position cannot overflow.
        let huge = format!("{{{ORDER}, \"position\": {{\"nu\": {}}}}}", u64::MAX);
        assert!(Header::parse(&huge).is_err());
    }

    #[test]
    fn extras_follow_the_fixed_part_even_when_it_holds_separators() {
        let h = Header::parse(&format!(
            "{{{ORDER}, \"extra\": [\"x\", \"deltap_air\", \"n_self\"], \"extra_separator\": \";\"}}"
        ))
        .unwrap();
        assert_eq!(
            h.optional,
            [Some(Source::Extra(2)), None, Some(Source::Extra(1))]
        );
        let fixed = format!("{:<160}", " 21  100.000000 1.000E-20;;;");
        let v = h.optional_values(&format!("{fixed};1;2e-5; 0.7 ;9"));
        assert_eq!(v[0], 0.7);
        assert!(v[1].is_nan());
        assert_eq!(v[2], 2e-5);
        let v = h.optional_values(&format!("{fixed};1;#"));
        assert!(v.iter().all(|x| x.is_nan()));
        assert!(h.optional_values(&fixed).iter().all(|x| x.is_nan()));
    }

    #[test]
    fn optional_fields_in_order_are_fixed_width_columns() {
        // As hapi.cache2storage writes them: after the 160 characters, without `position`.
        let order = ORDER
            .replace(r#""quanta"],"#, r#""quanta", "n_self", "deltap_air"],"#)
            .replace(
                r#""quanta": "%93s""#,
                r#""quanta": "%93s", "n_self": "%7.4f", "deltap_air": "%10.3E""#,
            );
        let h = Header::parse(&format!(
            "{{{order}, \"extra\": [\"n_self\", \"delta_self\"]}}"
        ))
        .unwrap();
        assert_eq!(
            h.optional,
            [
                Some(Source::Fixed(160..167)),
                Some(Source::Extra(1)),
                Some(Source::Fixed(167..177)),
            ]
        );
        assert_eq!(h.fixed_end, 177);
        let fixed = format!("{:<160}", " 21  100.000000 1.000E-20");
        let v = h.optional_values(&format!("{fixed} 0.6500 1.500D-05,#,-0.002"));
        assert_eq!(v, [0.65, -0.002, 1.5e-5]);
        let v = h.optional_values(&format!("{fixed}    #  garbage!!"));
        assert!(v.iter().all(|x| x.is_nan()), "{v:?}");
        // A record too short for a column (here n_self, 7 wide) lacks the value.
        let v = h.optional_values(&format!("{fixed} 0.65"));
        assert!(v.iter().all(|x| x.is_nan()), "{v:?}");
    }
}
