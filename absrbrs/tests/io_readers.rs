//! The line-file readers, on synthetic records written here (HITRAN data is never committed).

use std::fmt::Write as _;
use std::path::{Path, PathBuf};

use absrbrs::Error;
use absrbrs::io::{read_hapi, read_hitran_par};
use absrbrs::lines::{LineField, LineList};

/// `v` with `prec` decimals right-aligned in `width` characters, dropping the leading zero when
/// needed, as HITRAN's Fortran formats do (".0712", "-.004500").
fn fixed(v: f64, width: usize, prec: usize) -> String {
    let mut s = format!("{v:.prec$}");
    if s.len() > width {
        s = s.replacen("0.", ".", 1);
    }
    assert!(s.len() <= width, "{s} does not fit in {width}");
    format!("{s:>width$}")
}

/// One synthetic 160-character HITRAN record (the quanta, error codes, references and
/// statistical weights, bytes 67..160, are blank: the readers do not use them).
#[allow(clippy::too_many_arguments)] // one argument per field of the record
fn record(
    mol: u16,
    iso: char,
    nu: f64,
    s: f64,
    gair: f64,
    gself: f64,
    el: f64,
    nair: f64,
    dair: f64,
) -> String {
    let mut r = String::new();
    write!(
        r,
        "{mol:2}{iso}{}{s:10.3E}{:10.3E}{}{}{}{}{}",
        fixed(nu, 12, 6),
        1.0,
        fixed(gair, 5, 4),
        fixed(gself, 5, 3),
        fixed(el, 10, 4),
        fixed(nair, 4, 2),
        fixed(dair, 8, 6)
    )
    .unwrap();
    assert_eq!(r.len(), 67, "{r:?}");
    format!("{r:<160}")
}

/// A temporary directory, removed when dropped.
struct TempDir(PathBuf);

impl std::ops::Deref for TempDir {
    type Target = Path;
    fn deref(&self) -> &Path {
        &self.0
    }
}

impl AsRef<Path> for TempDir {
    fn as_ref(&self) -> &Path {
        &self.0
    }
}

impl Drop for TempDir {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// A fresh temporary directory; `name` must be unique among the tests of this file.
fn tmp(name: &str) -> TempDir {
    let dir = std::env::temp_dir().join(format!("absrbrs-io-{}-{name}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    TempDir(dir)
}

/// Writes `text` as `<dir>/f.par` and reads it.
fn read_par_text(dir: &Path, text: impl AsRef<[u8]>) -> absrbrs::Result<LineList> {
    let path = dir.join("f.par");
    std::fs::write(&path, text).unwrap();
    read_hitran_par(&path)
}

/// The line number and message of a parse error, or a panic for anything else.
fn parse_error(result: absrbrs::Result<LineList>) -> (usize, String) {
    match result {
        Err(Error::Parse {
            line_number,
            message,
            ..
        }) => (line_number, message),
        other => panic!("expected a parse error, got {other:?}"),
    }
}

/// A valid record with ν₀ = 6300.5.
fn good() -> String {
    record(2, '1', 6300.5, 1e-23, 0.07, 0.09, 1.0, 0.75, 0.0)
}

#[test]
fn reads_par_records_and_isotopologue_codes() {
    let dir = tmp("par");
    let path = dir.join("CO2.par");
    let text = [
        record(
            2, '1', 6300.5, 1.234e-23, 0.0712, 0.093, 123.4567, 0.75, -0.0045,
        ),
        record(2, '0', 6200.25, 2.5e-25, 0.07, 0.09, 50.0, 0.7, -0.004),
        record(2, 'A', 6250.0, 3.0e-26, 0.07, 0.09, 50.0, 0.7, 0.0),
        record(2, '2', 6100.0, 0.0, 0.07, 0.09, 50.0, 0.7, 0.0), // S = 0: skipped
    ]
    .join("\n");
    std::fs::write(&path, text + "\n").unwrap();
    let lines = read_hitran_par(&path).unwrap();
    assert_eq!(lines.len(), 3);
    assert_eq!(
        lines.column(LineField::Nu).unwrap(),
        &[6200.25, 6250.0, 6300.5]
    );
    assert_eq!(lines.isotopologue(), &[10, 11, 1]);
    assert_eq!(lines.molecule(), &[2, 2, 2]);
    assert_eq!(lines.value(LineField::SRef, 2), 1.234e-23);
    assert_eq!(lines.value(LineField::GammaAir, 2), 0.0712);
    assert_eq!(lines.value(LineField::ELower, 2), 123.4567);
    assert_eq!(lines.value(LineField::DeltaAir, 2), -0.0045);
    assert!(lines.column(LineField::NSelf).is_none());
}

#[test]
fn malformed_par_records_are_parse_errors_with_line_numbers() {
    let dir = tmp("bad");
    let path = dir.join("bad.par");
    let good = record(2, '1', 6300.5, 1e-23, 0.07, 0.09, 1.0, 0.75, 0.0);
    let mut bad = good.clone();
    bad.replace_range(3..15, "   abc.defgh");
    std::fs::write(&path, format!("{good}\n{bad}\n")).unwrap();
    let err = read_hitran_par(&path).unwrap_err();
    assert!(
        matches!(err, absrbrs::Error::Parse { line_number: 2, .. }),
        "{err}"
    );
    std::fs::write(&path, format!("{good}\n{}\n", &good[..100])).unwrap();
    assert!(matches!(
        read_hitran_par(&path),
        Err(absrbrs::Error::Parse { line_number: 2, .. })
    ));
}

#[test]
fn missing_file_is_an_io_error() {
    let err = read_hitran_par("/nonexistent/absrbrs/CO2.par").unwrap_err();
    assert!(matches!(
        err,
        absrbrs::Error::Io {
            kind: std::io::ErrorKind::NotFound,
            ..
        }
    ));
}

/// A HAPI header with the default fixed part and extra columns, as `hapi.fetch` writes it.
fn hapi_header(extra: &[&str]) -> String {
    let order = [
        "molec_id",
        "local_iso_id",
        "nu",
        "sw",
        "a",
        "gamma_air",
        "gamma_self",
        "elower",
        "n_air",
        "delta_air",
    ];
    let format = [
        ("molec_id", "%2d"),
        ("local_iso_id", "%1d"),
        ("nu", "%12.6f"),
        ("sw", "%10.3E"),
        ("a", "%10.3E"),
        ("gamma_air", "%5.4f"),
        ("gamma_self", "%5.3f"),
        ("elower", "%10.4f"),
        ("n_air", "%4.2f"),
        ("delta_air", "%8.6f"),
    ];
    let position = [
        ("molec_id", 0),
        ("local_iso_id", 2),
        ("nu", 3),
        ("sw", 15),
        ("a", 25),
        ("gamma_air", 35),
        ("gamma_self", 40),
        ("elower", 45),
        ("n_air", 55),
        ("delta_air", 59),
    ];
    let extra_format: Vec<String> = extra
        .iter()
        .map(|e| format!("\"{e}\": \"%9.6f\""))
        .collect();
    format!(
        "{{\"table_name\": \"T\", \"order\": {order:?}, \"format\": {{{}}}, \"position\": {{{}}}, \"extra\": {extra:?}, \"extra_format\": {{{}}}, \"extra_separator\": \",\"}}",
        format
            .iter()
            .map(|(k, v)| format!("\"{k}\": \"{v}\""))
            .collect::<Vec<_>>()
            .join(", "),
        position
            .iter()
            .map(|(k, v)| format!("\"{k}\": {v}"))
            .collect::<Vec<_>>()
            .join(", "),
        extra_format.join(", ")
    )
}

#[test]
fn reads_hapi_caches_with_optional_extras() {
    let dir = tmp("hapi");
    std::fs::write(
        dir.join("CO.header"),
        hapi_header(&["gamma_sdv_0_air_296", "n_self", "delta_self"]),
    )
    .unwrap();
    let a = record(5, '1', 2140.0, 1e-19, 0.05, 0.06, 10.0, 0.7, -0.002);
    let b = record(5, '2', 2100.0, 1e-21, 0.05, 0.06, 10.0, 0.7, -0.002);
    std::fs::write(
        dir.join("CO.data"),
        format!("{a},0.0500,0.65,-0.001\n{b},#,#,0.003\n"),
    )
    .unwrap();
    let lines = read_hapi(&dir, "CO").unwrap();
    assert_eq!(lines.column(LineField::Nu).unwrap(), &[2100.0, 2140.0]);
    assert_eq!(lines.value(LineField::NSelf, 1), 0.65);
    assert!(!lines.has_value(LineField::NSelf, 0)); // "#": absent, defaults to n_air
    assert_eq!(lines.value(LineField::NSelf, 0), 0.7);
    assert_eq!(lines.value(LineField::DeltaSelf, 0), 0.003);
    assert!(lines.column(LineField::DeltapAir).is_none());
}

#[test]
fn hapi_without_extras_and_bad_headers() {
    let dir = tmp("hapi2");
    std::fs::write(dir.join("X.header"), hapi_header(&[])).unwrap();
    std::fs::write(
        dir.join("X.data"),
        record(7, '1', 13100.0, 1e-23, 0.04, 0.04, 1.0, 0.7, 0.0) + "\n",
    )
    .unwrap();
    assert_eq!(read_hapi(&dir, "X").unwrap().len(), 1);
    std::fs::write(dir.join("Y.header"), "{not json").unwrap();
    std::fs::write(dir.join("Y.data"), "").unwrap();
    assert!(matches!(
        read_hapi(&dir, "Y"),
        Err(absrbrs::Error::Parse { .. })
    ));
    assert!(matches!(
        read_hapi(&dir, "missing"),
        Err(absrbrs::Error::Io { .. })
    ));
}

#[test]
fn crlf_endings_and_blank_lines_are_accepted() {
    let dir = tmp("crlf");
    let a = good();
    let b = record(2, '1', 6301.0, 1e-23, 0.07, 0.09, 1.0, 0.75, 0.0);
    // CRLF endings, an empty line, a whitespace-only line, and a final lone CR.
    let text = format!("{a}\r\n\r\n\n   \t\r\n{b}\r\n\r");
    assert_eq!(read_par_text(&dir, text).unwrap().len(), 2);
    let text = format!("\n\n{a}\n\n\n{b}");
    assert_eq!(read_par_text(&dir, text).unwrap().len(), 2);
}

#[test]
fn cr_only_endings_are_parse_errors() {
    let dir = tmp("cr");
    let (line, message) = parse_error(read_par_text(&dir, format!("{}\r{}\r", good(), good())));
    assert_eq!(line, 1);
    assert!(message.contains("carriage return"), "{message}");
    let text = format!("{}\n{}\r{}\n", good(), good(), good());
    assert_eq!(parse_error(read_par_text(&dir, text)).0, 2);
}

#[test]
fn non_ascii_records_are_parse_errors_with_line_numbers() {
    let dir = tmp("ascii");
    let mut bad = good();
    bad.replace_range(150..152, "é");
    let (line, message) = parse_error(read_par_text(&dir, format!("{}\n\n{bad}\n", good())));
    assert_eq!(line, 3);
    assert!(message.contains("non-ASCII"), "{message}");
    // Not UTF-8 either.
    let mut bytes = format!("{}\n", good()).into_bytes();
    bytes.extend_from_slice(&[0xff; 160]);
    assert_eq!(parse_error(read_par_text(&dir, bytes)).0, 2);
}

#[test]
fn negative_intensities_and_widths_are_parse_errors_with_line_numbers() {
    let dir = tmp("negative");
    for (range, text, field) in [
        (15..25, "-1.000E-23", "s_ref"),
        (35..40, "-.010", "gamma_air"),
        (40..45, "-0.01", "gamma_self"),
    ] {
        let mut bad = good();
        bad.replace_range(range, text);
        let (line, message) = parse_error(read_par_text(&dir, format!("{}\n{bad}\n", good())));
        assert_eq!(line, 2);
        assert!(
            message.contains(&format!("{field} is negative")),
            "{message}"
        );
    }
}

#[test]
fn fortran_number_forms_are_accepted() {
    let dir = tmp("fortran");
    let mut r = good();
    r.replace_range(15..25, "  1.234-23"); // a sign-only exponent
    r.replace_range(45..55, "  1.2345D2"); // a D exponent
    r.replace_range(59..67, "-4.50d-3");
    let lines = read_par_text(&dir, r).unwrap();
    assert_eq!(lines.value(LineField::SRef, 0), 1.234e-23);
    assert_eq!(lines.value(LineField::ELower, 0), 123.45);
    assert_eq!(lines.value(LineField::DeltaAir, 0), -4.5e-3);
    for not_finite in ["       nan", "       inf", "   1.0D400"] {
        let mut r = good();
        r.replace_range(15..25, not_finite);
        let (_, message) = parse_error(read_par_text(&dir, r));
        assert!(message.starts_with("s_ref: cannot parse"), "{message}");
    }
}

#[test]
fn hapi_missing_data_is_an_io_error() {
    let dir = tmp("nodata");
    std::fs::write(dir.join("X.header"), hapi_header(&[])).unwrap();
    let err = read_hapi(&dir, "X").unwrap_err();
    match err {
        Error::Io { path, kind, .. } => {
            assert_eq!(kind, std::io::ErrorKind::NotFound);
            assert!(path.ends_with("X.data"), "{path}");
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn hapi_extra_separator_defaults_to_comma() {
    let dir = tmp("separator");
    let header = hapi_header(&["n_self"]).replace(", \"extra_separator\": \",\"", "");
    assert!(!header.contains("extra_separator"));
    std::fs::write(dir.join("X.header"), header).unwrap();
    std::fs::write(dir.join("X.data"), format!("{},0.61\n", good())).unwrap();
    assert_eq!(
        read_hapi(&dir, "X").unwrap().value(LineField::NSelf, 0),
        0.61
    );
}

#[test]
fn hapi_non_utf8_header_is_a_parse_error() {
    let dir = tmp("utf8");
    std::fs::write(dir.join("X.header"), b"{\n\"order\": \xff}").unwrap();
    std::fs::write(dir.join("X.data"), "").unwrap();
    let (line, message) = parse_error(read_hapi(&dir, "X"));
    assert_eq!(line, 2);
    assert!(message.contains("UTF-8"), "{message}");
}

#[test]
fn hapi_optional_fields_in_order_are_fixed_width_columns() {
    // As hapi.cache2storage writes a table: the optional fields are fixed-width columns after
    // the 160 characters, and there is no `position`.
    let dir = tmp("storage");
    let header = serde_json::json!({
        "table_name": "T",
        "order": ["molec_id", "local_iso_id", "nu", "sw", "a", "gamma_air", "gamma_self",
            "elower", "n_air", "delta_air", "quanta", "n_self", "delta_self", "deltap_air"],
        "format": {"molec_id": "%2d", "local_iso_id": "%1d", "nu": "%12.6f", "sw": "%10.3E",
            "a": "%10.3E", "gamma_air": "%5.4f", "gamma_self": "%5.3f", "elower": "%10.4f",
            "n_air": "%4.2f", "delta_air": "%8.6f", "quanta": "%93s", "n_self": "%7.4f",
            "delta_self": "%9.6f", "deltap_air": "%10.3E"},
    });
    std::fs::write(dir.join("S.header"), header.to_string()).unwrap();
    let b = record(2, '1', 6301.0, 1e-23, 0.07, 0.09, 1.0, 0.75, 0.0);
    let data = format!(
        "{} 0.6100-0.001500 1.500E-05\n{b}      #        #          #\n",
        good()
    );
    std::fs::write(dir.join("S.data"), data).unwrap();
    let lines = read_hapi(&dir, "S").unwrap();
    assert_eq!(lines.column(LineField::NSelf).unwrap()[0], 0.61);
    assert_eq!(lines.column(LineField::DeltaSelf).unwrap()[0], -0.0015);
    assert_eq!(lines.column(LineField::DeltapAir).unwrap()[0], 1.5e-5);
    for field in [LineField::NSelf, LineField::DeltaSelf, LineField::DeltapAir] {
        assert!(!lines.has_value(field, 1), "{field:?}");
    }
}

/// Whether a reader's result is one it may return for arbitrary input: a line list or a
/// `Parse` error (the readers check every value per record, so `from_columns`'s
/// `InvalidInput` is unreachable).
fn acceptable(result: &absrbrs::Result<LineList>) -> bool {
    matches!(result, Ok(_) | Err(Error::Parse { .. }))
}

proptest::proptest! {
    #![proptest_config(proptest::prelude::ProptestConfig::with_cases(512))]
    #[test]
    fn par_reader_never_panics(text in "(?s).{0,300}") {
        let dir = tmp("fuzz");
        let result = read_par_text(&dir, &text);
        proptest::prop_assert!(acceptable(&result), "unexpected {result:?}");
    }
}

proptest::proptest! {
    // Each case writes three files, so fewer cases keep the test quick.
    #![proptest_config(proptest::prelude::ProptestConfig::with_cases(256))]
    /// A valid record with bytes overwritten in and around its fields, then truncated or
    /// extended, read as a `.par` record and as a HAPI record with an extra column.
    #[test]
    fn mutated_records_never_panic(
        edits in proptest::collection::vec((0usize..80, 0u8..0x80), 0..8),
        len in 0usize..200,
        fill in 0x20u8..0x7f,
    ) {
        let dir = tmp("mutated");
        let mut bytes = good().into_bytes();
        bytes.resize(len, fill);
        for (at, byte) in edits {
            if let Some(b) = bytes.get_mut(at) {
                *b = byte;
            }
        }
        let result = read_par_text(&dir, &bytes);
        proptest::prop_assert!(acceptable(&result), "unexpected {result:?}");
        std::fs::write(dir.join("H.header"), hapi_header(&["n_self"])).unwrap();
        bytes.extend_from_slice(b",0.5\n");
        std::fs::write(dir.join("H.data"), &bytes).unwrap();
        let result = read_hapi(&dir, "H");
        proptest::prop_assert!(acceptable(&result), "unexpected {result:?}");
    }
}
