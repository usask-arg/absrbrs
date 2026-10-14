//! HAPI caches: a JSON `.header` describing a `.data` file of fixed-width records followed by
//! separated extra columns.

use std::collections::HashMap;
use std::ops::Range;
use std::path::Path;

use serde_json::{Map, Value};

use super::hitran::{FixedLayout, parse_fixed};
use super::{Builder, Lines, io_error, parse_error};
use crate::Result;
use crate::lines::LineList;

/// The optional fields read from the extras, in [`Builder`] order.
const OPTIONAL: [&str; 3] = ["n_self", "delta_self", "deltap_air"];

/// What the reader needs from a `.header`.
#[derive(Debug)]
struct Header {
    /// Byte ranges of the fixed fields that are read.
    layout: FixedLayout,
    /// Where the fixed-width part of a record ends (the extras follow the next separator).
    fixed_end: usize,
    separator: String,
    /// The index in `extra` of each of [`OPTIONAL`].
    extra_index: [Option<usize>; 3],
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
        let extra_index = OPTIONAL.map(|name| extra.iter().position(|&e| e == name));
        Ok(Header {
            layout,
            fixed_end,
            separator,
            extra_index,
        })
    }

    /// The n_self, δ_self and δ′_air of a record: NaN when absent from the header, missing from
    /// the record, or not a finite number (HAPI writes `#` for missing).
    fn optional_values(&self, record: &str) -> [f64; 3] {
        let mut values = [f64::NAN; 3];
        let Some(last) = self.extra_index.iter().flatten().max().copied() else {
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
            for (value, index) in values.iter_mut().zip(self.extra_index) {
                if index == Some(k) {
                    *value = token
                        .trim()
                        .parse::<f64>()
                        .ok()
                        .filter(|v| v.is_finite())
                        .unwrap_or(f64::NAN);
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
/// δ′_air are read from the extra columns when the header's `extra` lists them; a missing value
/// (`#`) is NaN, meaning absent for that line. Without them in `extra`, the column is absent.
///
/// Records with S_ref = 0 or ν₀ ≤ 0 are skipped, and so are empty lines.
///
/// # Errors
///
/// - [`Error::Io`](crate::Error::Io) if either file cannot be read;
/// - [`Error::Parse`](crate::Error::Parse) for a header that is not JSON or lacks a field
///   (line 1 of the header), and, with the line number, for a data record that is not ASCII,
///   is too short for its fixed fields, or has a fixed field that does not parse;
/// - [`Error::InvalidInput`](crate::Error::InvalidInput) if the values do not make a valid
///   [`LineList`], for example a negative width.
pub fn read_hapi(directory: impl AsRef<Path>, table: &str) -> Result<LineList> {
    let directory = directory.as_ref();
    let header_path = directory.join(format!("{table}.header"));
    let shown = header_path.display().to_string();
    let text = std::fs::read_to_string(&header_path).map_err(|e| io_error(&shown, &e))?;
    let header = Header::parse(&text).map_err(|message| parse_error(&shown, 1, message))?;

    let mut lines = Lines::open(&directory.join(format!("{table}.data")))?;
    let path = lines.path.clone();
    let min_len = header.layout.min_len();
    let mut builder = Builder::new(header.extra_index.map(|k| k.is_some()));
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
        if fields.is_kept() {
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
        assert_eq!(h.extra_index, [None; 3]);
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
        assert_eq!(h.extra_index, [Some(2), None, Some(1)]);
        let fixed = format!("{:<160}", " 21  100.000000 1.000E-20;;;");
        let v = h.optional_values(&format!("{fixed};1;2e-5; 0.7 ;9"));
        assert_eq!(v[0], 0.7);
        assert!(v[1].is_nan());
        assert_eq!(v[2], 2e-5);
        let v = h.optional_values(&format!("{fixed};1;#"));
        assert!(v.iter().all(|x| x.is_nan()));
        assert!(h.optional_values(&fixed).iter().all(|x| x.is_nan()));
    }
}
