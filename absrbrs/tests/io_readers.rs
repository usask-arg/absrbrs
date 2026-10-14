//! The line-file readers, on synthetic records written here (HITRAN data is never committed).

use std::fmt::Write as _;
use std::path::PathBuf;

use absrbrs::io::{read_hapi, read_hitran_par};
use absrbrs::lines::LineField;

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

fn tmp(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("absrbrs-io-{}-{name}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    dir
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

proptest::proptest! {
    #![proptest_config(proptest::prelude::ProptestConfig::with_cases(512))]
    #[test]
    fn par_reader_never_panics(text in "(?s).{0,300}") {
        let dir = tmp("fuzz");
        let path = dir.join("f.par");
        std::fs::write(&path, &text).unwrap();
        match read_hitran_par(&path) {
            Ok(_) | Err(absrbrs::Error::Parse { .. }) | Err(absrbrs::Error::InvalidInput(_)) => {}
            Err(e) => proptest::prop_assert!(false, "unexpected error {e}"),
        }
    }
}
