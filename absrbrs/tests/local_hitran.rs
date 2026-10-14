//! Opt-in: read the user's own HAPI caches (e.g. sasktran2's database/hitran_lines). Set
//! ABSRBRS_HITRAN_LINES_DIR to that directory; skipped otherwise (and never run in CI, since
//! HITRAN data is not redistributed).

use std::io::BufRead;
use std::path::Path;

use absrbrs::io::read_hapi;
use absrbrs::lines::{LineField, isotopologue};

/// Counts the records of `<table>.data`: (skipped for ν₀ ≤ 0, skipped for S_ref = 0, all).
/// Reads ν₀ and S_ref at their default HITRAN offsets (bytes 3..15 and 15..25), which the
/// headers of these caches use.
fn skipped_records(dir: &Path, table: &str) -> (usize, usize, usize) {
    let file = std::fs::File::open(dir.join(format!("{table}.data"))).unwrap();
    let (mut nu_le_0, mut s_zero, mut records) = (0, 0, 0);
    for line in std::io::BufReader::new(file).split(b'\n') {
        let line = line.unwrap();
        if line.is_empty() || line == b"\r" {
            continue;
        }
        records += 1;
        let field = |r: std::ops::Range<usize>| -> f64 {
            std::str::from_utf8(&line[r])
                .unwrap()
                .trim()
                .parse()
                .unwrap()
        };
        if field(15..25) == 0.0 {
            s_zero += 1;
        } else if field(3..15) <= 0.0 {
            nu_le_0 += 1;
        }
    }
    (nu_le_0, s_zero, records)
}

#[test]
fn reads_local_hapi_caches() {
    let Ok(dir) = std::env::var("ABSRBRS_HITRAN_LINES_DIR") else {
        eprintln!("skipped: set ABSRBRS_HITRAN_LINES_DIR to a directory of HAPI caches");
        return;
    };
    let mut tables: Vec<_> = std::fs::read_dir(&dir)
        .unwrap()
        .filter_map(|e| {
            e.ok()?
                .path()
                .file_name()?
                .to_str()?
                .strip_suffix(".header")
                .map(String::from)
        })
        .collect();
    tables.sort();
    for table in tables {
        let t = std::time::Instant::now();
        let lines = read_hapi(&dir, &table).unwrap();
        let nu = lines.column(LineField::Nu).unwrap();
        assert!(nu.windows(2).all(|w| w[0] <= w[1]));
        let elapsed = t.elapsed();
        let unknown = (0..lines.len())
            .filter(|&i| isotopologue(lines.molecule()[i], lines.isotopologue()[i]).is_none())
            .count();
        let (nu_le_0, s_zero, records) = skipped_records(Path::new(&dir), &table);
        assert_eq!(records - nu_le_0 - s_zero, lines.len());
        eprintln!(
            "{table}: {} lines in {elapsed:.2?}; skipped {nu_le_0} with nu <= 0 and {s_zero} \
             with S = 0; {unknown} without isotopologue data",
            lines.len(),
        );
    }
}
