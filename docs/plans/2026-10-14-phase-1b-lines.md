# absrbrs Phase 1b (Line Data, Partition Functions, Line State) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use usask-foundation:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the inputs of a line-by-line cross section:
- line lists, read from HITRAN `.par` and HAPI caches;
- isotopologue masses and abundances;
- built-in TIPS-2021 partition functions, with an analytic d ln Q/dT;
- the line-state layer, which maps a condition and a line to q = (A, ν₀′, σ_D, y) with closed-form ∂q/∂θ.

Python gets `absrbrs.LineList` and `absrbrs.partition`. This is phase 1b of `docs/specs/2026-10-14-absrbrs-design.md` (§5, §13).

**Architecture:** Five new pieces.
- **`absrbrs::lines`** holds a structure-of-arrays `LineList`, always sorted by ν₀, and the isotopologue table. The table is generated from HAPI's MIT-licensed `hapi.ISO`.
- **`absrbrs::io`** has the readers.
- **`absrbrs::partition`** has:
  - a `PartitionFunction` trait;
  - `Tips2021`, built from all 181 TIPS-2021 tables at 1 K steps, stored as f32 in an embedded file of about 4 MB;
  - `TabulatedPartition` for user tables.

  Both interpolate ln Q with cubic Hermite polynomials, so d ln Q/dT is the interpolant's own derivative.
- **`absrbrs::state`** (crate-private) is generic over a small `Real` trait:
  - it computes q in f64 for production;
  - two separate functions return the closed-form Jacobians, one for the state (T, p, VMR) and one for the ten line fields;
  - a test-only dual number checks them;
  - an Arb reference checks the values.
- **Python bindings** for line lists and partition functions.

**Tech Stack:**
- Rust (stable, edition 2024, MSRV 1.85) with num-complex and thiserror; new: serde_json (core).
- pyo3 and numpy 0.29.
- python-flint 0.9.0 (Arb); hitran-api 1.3.0.0, used only by the generator script and validation.
- uv inline-script tools.

## Global Constraints

These come from the spec, `docs/plans/phase-0-followups.md`, `docs/plans/phase-1a-followups.md` and the phase 1b design. Every task implicitly includes them.

- **Core API rules:**
  - no pyo3 or ndarray in the core public API;
  - `#![forbid(unsafe_code)]` and `#![warn(missing_docs)]`;
  - no panics on finite input;
  - one `#[non_exhaustive]` `Error` enum. A new variant needs mappings in `absrbrs-py/src/errors.rs`, `python/absrbrs/errors.py` and the `_testing` hook.
- **Units in the core API:**
  - wavenumber cm⁻¹, temperature K, pressure Pa, VMR dimensionless;
  - line parameters in HITRAN units: S_ref in cm⁻¹/(molecule cm⁻²); widths and shifts in cm⁻¹ atm⁻¹; δ′_air in cm⁻¹ atm⁻¹ K⁻¹;
  - T_ref = 296 K, 1 atm = 101 325 Pa;
  - masses in g/mol;
  - A in m²/molecule, so that σ(ν) = A·K(x, y).
- **Physical constants** are the exact SI values:
  - c = 299 792 458 m/s;
  - k_B = 1.380649 × 10⁻²³ J/K;
  - N_A = 6.02214076 × 10²³ mol⁻¹;
  - h = 6.62607015 × 10⁻³⁴ J s;
  - c₂ = hc/k_B = 1.4387768775039338 cm K (the nearest double; sasktran2 uses 1.4387769).
- **Isotopologue numbering** follows HITRAN. In `.par` records the local isotopologue character `1`–`9` is 1–9, `0` is 10, and `A`, `B`, … are 11, 12, ….
- **Optional line fields:** n_self (defaults to n_air), δ_self (defaults to δ_air) and δ′_air (defaults to 0).
  - A NaN in an optional column means "absent for this line".
  - **Derivatives with respect to an absent field (refines spec §5.3).** The derivative with respect to the absent field itself is the partial, taken at its default value. The derivative with respect to the field it defaults to *includes* its effect through the default. So with n_self absent, ∂/∂n_air contains both the air and the self terms, which keeps it consistent with finite differences of `with_column("n_air", …)`.
- **Licences** (`THIRD_PARTY.md`, rows confirmed 2026-10-14):
  - TIPS-2021 is MIT, copyright Robert Gamache. Build it offline from the Zenodo `.QTpy` files. Check their MD5s and use a restricted unpickler that refuses all globals; never unpickle at runtime. Ship `absrbrs/LICENSES/TIPS-2021-MIT.txt` and list it in `pyproject.toml` `license-files`.
  - HAPI is MIT, copyright 2018 HITRAN team. The isotopologue table derived from `hapi.ISO` is embedded, so ship `absrbrs/LICENSES/HAPI-MIT.txt` and list it in `license-files`.
  - HITRAN line data is **blocked**: never commit, ship or upload HITRAN records, including in CI artifacts. Tests use synthetic records written by the tests themselves. An opt-in local test may read the user's own caches.
- **Dependencies:** every new one needs the usask-foundation:dependency-audit skill. Phase 1b adds `serde_json` (core, moved up from dev).
- **Copying from scatrs:** only with `git -C ~/dev/scat-rs show 4793160:<path>`; never touch its working tree.
- **Formatting:** run `cargo fmt --all` and `uv run --locked ruff format .` before checks.
- **Clippy on current stable:** lint with both `cargo clippy` and `cargo +1.99.0 clippy --workspace --all-targets -- -D warnings`. CI uses the latest stable, and the local default is older.
- **Commits:** stage files explicitly (no `git commit -a`). Each ends with a `Co-Authored-By:` trailer naming the model that wrote it. Do not push.

## File Map

```
tools/isotopologues.py                       generator (PEP 723 script)              Task 1
absrbrs/src/lines/{mod.rs, isotopologues.rs, isotopologue_table.rs}                  Tasks 1, 2
absrbrs/LICENSES/HAPI-MIT.txt                                                        Task 1
absrbrs/src/io/{mod.rs, hitran.rs, hapi.rs}                                          Task 3
absrbrs/tests/{io_readers.rs, local_hitran.rs}                                       Task 3
tools/tips2021.py                            converter (PEP 723 script)              Task 4
absrbrs/data/{tips2021.bin, tips2021.json}   embedded tables + provenance            Task 4
absrbrs/src/partition/{mod.rs, hermite.rs, tips2021.rs, tabulated.rs}                Task 4
absrbrs/LICENSES/TIPS-2021-MIT.txt, tests/test_data_files.py                         Task 4
absrbrs/src/state/{mod.rs, real.rs, jacobian.rs, dual.rs}                            Task 5
reference/absrbrs_reference/generators/line_state.py, reference/data/line_state.jsonl, absrbrs/tests/line_state_reference.rs   Task 5
absrbrs-py/src/{lines.rs, partition.rs}, python/absrbrs/{lines.py, partition.py}, stubs, tests/test_lines.py, tests/test_partition.py   Task 6
validation/tests/test_tips.py, docs/theory/{line-state.md, partition-functions.md}, docs/plans/phase-1b-followups.md   Task 7
```

Work on branch `phase-1b-lines`, which already exists and is checked out.

---

### Task 1: Isotopologue table

**Files:**
- Create: `tools/isotopologues.py`, `absrbrs/src/lines/mod.rs`, `absrbrs/src/lines/isotopologues.rs`, `absrbrs/src/lines/isotopologue_table.rs` (generated), `absrbrs/LICENSES/HAPI-MIT.txt`
- Modify: `absrbrs/src/lib.rs` (`pub mod lines;`), `pyproject.toml` (`license-files`), `THIRD_PARTY.md` (HAPI row)

**Interfaces:**
- Produces, in `absrbrs::lines`:
  - `pub struct Isotopologue { pub molecule: u16, pub isotopologue: u16, pub global_id: u16, pub formula: &'static str, pub molecule_name: &'static str, pub abundance: f64, pub mass: f64 }`, with mass in g/mol;
  - `pub fn isotopologue(molecule: u16, isotopologue: u16) -> Option<&'static Isotopologue>`;
  - `pub fn molecule_id(name: &str) -> Option<u16>`;
  - `pub fn molecule_name(molecule: u16) -> Option<&'static str>`;
  - `pub fn isotopologues() -> &'static [Isotopologue]`.

- [ ] **Step 1: Write the generator**

Create `tools/isotopologues.py`:
```python
# /// script
# requires-python = ">=3.11"
# dependencies = ["hitran-api==1.3.0.0", "numpy"]
# ///
"""Generate absrbrs/src/lines/isotopologue_table.rs from HAPI's isotopologue table.

HAPI (hitran-api 1.3.0.0, MIT licence, Copyright 2018 HITRAN team) holds, for each HITRAN
(molecule, isotopologue) pair, its global id, formula, natural abundance and mass (g/mol) in
`hapi.ISO`. Run with `uv run tools/isotopologues.py`; the output is committed.
"""

from __future__ import annotations

import contextlib
import io
import pathlib
import re

with contextlib.redirect_stdout(io.StringIO()):  # HAPI prints a banner on import
    import hapi

OUT = pathlib.Path(__file__).resolve().parents[1] / "absrbrs" / "src" / "lines" / "isotopologue_table.rs"


def rust_float(v: float) -> str:
    """The shortest round-trip literal, without exponent zero-padding (clippy rejects e-08)."""
    s = repr(float(v))
    s = re.sub(r"e([+-])0*(\d)", r"e\1\2", s).replace("e+", "e")
    return s if any(c in s for c in ".e") else s + ".0"


def rust_str(s: str) -> str:
    return '"' + s.replace("\\", "\\\\").replace('"', '\\"') + '"'


def main() -> None:
    rows = sorted(hapi.ISO.items())
    out = [
        "// Generated by tools/isotopologues.py from HAPI 1.3.0.0's isotopologue table (hapi.ISO).",
        "// HAPI: MIT licence, Copyright 2018 HITRAN team; see absrbrs/LICENSES/HAPI-MIT.txt.",
        "// Do not edit by hand.",
        "",
        "use super::isotopologues::Isotopologue;",
        "",
        "#[rustfmt::skip]",
        f"pub(crate) const ISOTOPOLOGUES: [Isotopologue; {len(rows)}] = [",
    ]
    for (m, i), (gid, formula, abundance, mass, mol) in rows:
        out.append(
            f"    Isotopologue {{ molecule: {m}, isotopologue: {i}, global_id: {gid}, "
            f"formula: {rust_str(formula)}, molecule_name: {rust_str(mol)}, "
            f"abundance: {rust_float(abundance)}, mass: {rust_float(mass)} }},"
        )
    out.append("];")
    OUT.write_text("\n".join(out) + "\n", encoding="utf-8")
    print(f"{OUT}: {len(rows)} isotopologues")


if __name__ == "__main__":
    main()
```
Run: `uv run tools/isotopologues.py`. Expected: `…/isotopologue_table.rs: 156 isotopologues`.

- [ ] **Step 2: Write the failing tests and the lookup module**

Create `absrbrs/src/lines/isotopologues.rs`:
```rust
//! Isotopologue metadata in HITRAN's numbering: masses and natural abundances (from HAPI's
//! table, MIT licence; see absrbrs/LICENSES/HAPI-MIT.txt).

use super::isotopologue_table::ISOTOPOLOGUES;

/// One isotopologue: HITRAN molecule and (local) isotopologue numbers, and its properties.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Isotopologue {
    /// HITRAN molecule number (1 = H₂O, 2 = CO₂, …).
    pub molecule: u16,
    /// HITRAN local isotopologue number (1, 2, …; `.par` codes 0, A, B are 10, 11, 12).
    pub isotopologue: u16,
    /// HITRAN global isotopologue id.
    pub global_id: u16,
    /// Formula, e.g. `H2(16O)`.
    pub formula: &'static str,
    /// Molecule name, e.g. `H2O`.
    pub molecule_name: &'static str,
    /// Natural terrestrial abundance (included in HITRAN line intensities).
    pub abundance: f64,
    /// Molar mass, g/mol.
    pub mass: f64,
}

/// The isotopologue (molecule, isotopologue), if the table has it.
pub fn isotopologue(molecule: u16, isotopologue: u16) -> Option<&'static Isotopologue> {
    ISOTOPOLOGUES
        .binary_search_by_key(&(molecule, isotopologue), |i| (i.molecule, i.isotopologue))
        .ok()
        .map(|k| &ISOTOPOLOGUES[k])
}

/// The HITRAN molecule number of `name` (exact, case-sensitive: "CO2", "H2O", …).
pub fn molecule_id(name: &str) -> Option<u16> {
    ISOTOPOLOGUES.iter().find(|i| i.molecule_name == name).map(|i| i.molecule)
}

/// The name of HITRAN molecule `molecule`.
pub fn molecule_name(molecule: u16) -> Option<&'static str> {
    ISOTOPOLOGUES.iter().find(|i| i.molecule == molecule).map(|i| i.molecule_name)
}

/// Every isotopologue in the table, sorted by (molecule, isotopologue).
pub fn isotopologues() -> &'static [Isotopologue] {
    &ISOTOPOLOGUES
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn table_is_sorted_and_unique() {
        assert!(isotopologues().windows(2).all(|w| (w[0].molecule, w[0].isotopologue) < (w[1].molecule, w[1].isotopologue)));
        assert_eq!(isotopologues().len(), 156);
    }

    #[test]
    fn known_entries() {
        let h2o = isotopologue(1, 1).unwrap();
        assert_eq!((h2o.formula, h2o.molecule_name, h2o.abundance, h2o.mass), ("H2(16O)", "H2O", 0.9973173, 18.010565));
        let co2 = isotopologue(2, 12).unwrap();
        assert_eq!((co2.formula, co2.mass), ("(13C)(17O)2", 47.001618));
        assert!(isotopologue(2, 99).is_none());
    }

    #[test]
    fn molecule_names() {
        assert_eq!(molecule_id("CO2"), Some(2));
        assert_eq!(molecule_id("O2"), Some(7));
        assert_eq!(molecule_id("co2"), None);
        assert_eq!(molecule_name(1), Some("H2O"));
        assert_eq!(molecule_name(0), None);
    }
}
```
Create `absrbrs/src/lines/mod.rs`:
```rust
//! Spectroscopic line lists and isotopologue metadata (spec §5.1, §5.6).

mod isotopologue_table;
mod isotopologues;

pub use isotopologues::{Isotopologue, isotopologue, isotopologues, molecule_id, molecule_name};
```
Add `pub mod lines;` to `absrbrs/src/lib.rs`, after `pub mod faddeeva;`.

Run: `cargo test -p absrbrs lines`. Expected: `3 passed`. If `known_entries` fails on a value, the generator output differs from HAPI 1.3.0.0; check the pinned version. Do not edit the generated file.

- [ ] **Step 3: Licence notice and register**

1. Create `absrbrs/LICENSES/HAPI-MIT.txt`. It holds:
   - the heading `HAPI isotopologue table (absrbrs/src/lines/isotopologue_table.rs is derived from hapi.ISO)`;
   - a blank line;
   - HAPI's `LICENSE.txt` text, exactly as it appears in https://github.com/hitranonline/hapi (MIT, "Copyright 2018 HITRAN team"). Fetch it with `curl -s https://raw.githubusercontent.com/hitranonline/hapi/master/LICENSE.txt`, using no personal data in the request.
2. In `pyproject.toml`, set `license-files = ["LICENSE", "absrbrs/LICENSES/Faddeeva-MIT.txt", "absrbrs/LICENSES/HAPI-MIT.txt"]`.
3. In `THIRD_PARTY.md`'s HAPI row, change the use to "isotopologue table derived from `hapi.ISO`, embedded in the crate (`absrbrs/src/lines/isotopologue_table.rs`, generated by `tools/isotopologues.py`); optional runtime use for cross-code comparisons". Change the status to **conditions**, with the obligation "ship `absrbrs/LICENSES/HAPI-MIT.txt`, listed in `license-files`".

- [ ] **Step 4: Check and commit**

Run:
```bash
cargo test -p absrbrs
cargo clippy --workspace --all-targets -- -D warnings
cargo +1.99.0 clippy --workspace --all-targets -- -D warnings
cargo fmt --all -- --check
uv run --locked ruff check . && uv run --locked ruff format --check .
```
Expected: all clean. `ruff` also checks `tools/isotopologues.py`; run `ruff format` if needed.

Commit:
```bash
git add tools/isotopologues.py absrbrs/src/lib.rs absrbrs/src/lines absrbrs/LICENSES/HAPI-MIT.txt pyproject.toml THIRD_PARTY.md
git commit -m "feat(lines): isotopologue table generated from HAPI (masses, abundances)

Co-Authored-By: <your model> <noreply@anthropic.com>"
```

---

### Task 2: `LineList`

**Files:**
- Create: `absrbrs/src/lines/list.rs`
- Modify: `absrbrs/src/lines/mod.rs`

**Interfaces:**
- Consumes: `crate::{Error, Result}`.
- Produces, in `absrbrs::lines`:
  - `pub enum LineField { Nu, SRef, ELower, GammaAir, GammaSelf, NAir, NSelf, DeltaAir, DeltapAir, DeltaSelf }`
    - derives `Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord`;
    - `pub const ALL: [LineField; 10]`;
    - `pub fn name(self) -> &'static str` (`"nu"`, `"s_ref"`, `"e_lower"`, `"gamma_air"`, `"gamma_self"`, `"n_air"`, `"n_self"`, `"delta_air"`, `"deltap_air"`, `"delta_self"`);
    - `pub fn from_name(name: &str) -> Option<LineField>`;
    - `pub fn is_optional(self) -> bool`.
  - `pub struct LineColumns { pub molecule: Vec<u16>, pub isotopologue: Vec<u16>, pub nu, s_ref, e_lower, gamma_air, gamma_self, n_air, delta_air: Vec<f64>, pub n_self, delta_self, deltap_air: Option<Vec<f64>> }`, deriving `Clone, Debug, Default, PartialEq`.
  - `pub struct LineList`, deriving `Clone, Debug, PartialEq`, with:
    - `from_columns(LineColumns) -> Result<LineList>`
    - `len()`, `is_empty()`
    - `molecule() -> &[u16]`, `isotopologue() -> &[u16]`
    - `column(LineField) -> Option<&[f64]>` (None only for an absent optional column)
    - `value(LineField, usize) -> f64` (resolves defaults: n_self → n_air, δ_self → δ_air, δ′_air → 0)
    - `has_value(LineField, usize) -> bool`
    - `select(molecule: Option<u16>, isotopologues: Option<&[u16]>, wavenumber: Option<(f64, f64)>) -> LineList`
    - `with_column(LineField, Vec<f64>) -> Result<LineList>`
- Also produces a crate-private constructor, `LineList::from_unsorted(LineColumns)`, for use by the readers. It is the same as `from_columns`.

- [ ] **Step 1: Write the failing tests**

Create `absrbrs/src/lines/list.rs` with this test module at the end (the implementation goes above it in Step 3):
```rust
#[cfg(test)]
mod tests {
    use super::*;

    fn columns(nu: Vec<f64>) -> LineColumns {
        let n = nu.len();
        LineColumns {
            molecule: vec![2; n],
            isotopologue: (1..=n as u16).collect(),
            s_ref: (0..n).map(|i| 1e-21 * (i + 1) as f64).collect(),
            e_lower: vec![100.0; n],
            gamma_air: vec![0.07; n],
            gamma_self: vec![0.09; n],
            n_air: vec![0.75; n],
            delta_air: vec![-0.004; n],
            nu,
            ..Default::default()
        }
    }

    #[test]
    fn sorts_by_wavenumber_keeping_rows_together() {
        let lines = LineList::from_columns(columns(vec![3.0, 1.0, 2.0])).unwrap();
        assert_eq!(lines.column(LineField::Nu).unwrap(), &[1.0, 2.0, 3.0]);
        assert_eq!(lines.isotopologue(), &[2, 3, 1]);
        assert_eq!(lines.column(LineField::SRef).unwrap(), &[2e-21, 3e-21, 1e-21]);
        assert_eq!(lines.len(), 3);
    }

    #[test]
    fn optional_fields_default_per_line() {
        let mut c = columns(vec![1.0, 2.0]);
        c.n_self = Some(vec![0.6, f64::NAN]);
        let lines = LineList::from_columns(c).unwrap();
        assert_eq!(lines.value(LineField::NSelf, 0), 0.6);
        assert_eq!(lines.value(LineField::NSelf, 1), 0.75); // NaN: absent, so n_air
        assert!(lines.has_value(LineField::NSelf, 0) && !lines.has_value(LineField::NSelf, 1));
        assert_eq!(lines.value(LineField::DeltaSelf, 0), -0.004); // column absent: delta_air
        assert_eq!(lines.value(LineField::DeltapAir, 1), 0.0);
        assert!(lines.column(LineField::DeltaSelf).is_none());
    }

    #[test]
    fn rejects_bad_columns() {
        let mut c = columns(vec![1.0, 2.0]);
        c.gamma_air.pop();
        assert!(matches!(LineList::from_columns(c), Err(Error::InvalidInput(_))));
        let mut c = columns(vec![1.0, 2.0]);
        c.s_ref[1] = f64::NAN;
        assert!(matches!(LineList::from_columns(c), Err(Error::InvalidInput(_))));
        let mut c = columns(vec![1.0, 0.0]);
        c.nu[1] = 0.0;
        assert!(matches!(LineList::from_columns(c), Err(Error::InvalidInput(_))));
        let mut c = columns(vec![1.0, 2.0]);
        c.n_self = Some(vec![0.5, f64::INFINITY]);
        assert!(matches!(LineList::from_columns(c), Err(Error::InvalidInput(_))));
    }

    #[test]
    fn select_filters() {
        let mut c = columns(vec![1.0, 2.0, 3.0, 4.0]);
        c.molecule = vec![1, 2, 2, 2];
        let lines = LineList::from_columns(c).unwrap();
        let s = lines.select(Some(2), Some(&[3, 4]), Some((2.5, 4.0)));
        assert_eq!(s.column(LineField::Nu).unwrap(), &[3.0, 4.0]);
        assert_eq!(lines.select(None, None, None), lines);
        assert!(lines.select(Some(9), None, None).is_empty());
    }

    #[test]
    fn with_column_replaces_and_resorts() {
        let lines = LineList::from_columns(columns(vec![1.0, 2.0])).unwrap();
        let g = lines.with_column(LineField::GammaAir, vec![0.01, 0.02]).unwrap();
        assert_eq!(g.column(LineField::GammaAir).unwrap(), &[0.01, 0.02]);
        let r = lines.with_column(LineField::Nu, vec![5.0, 4.0]).unwrap();
        assert_eq!(r.column(LineField::Nu).unwrap(), &[4.0, 5.0]);
        assert_eq!(r.isotopologue(), &[2, 1]);
        let s = lines.with_column(LineField::NSelf, vec![0.5, f64::NAN]).unwrap();
        assert_eq!(s.value(LineField::NSelf, 1), 0.75);
        assert!(lines.with_column(LineField::NAir, vec![1.0]).is_err());
    }

    #[test]
    fn field_names_round_trip() {
        for f in LineField::ALL {
            assert_eq!(LineField::from_name(f.name()), Some(f));
        }
        assert_eq!(LineField::from_name("sw"), None);
        assert_eq!(LineField::ALL.iter().filter(|f| f.is_optional()).count(), 3);
    }
}
```
Add `mod list;` and `pub use list::{LineColumns, LineField, LineList};` to `absrbrs/src/lines/mod.rs`.

Run: `cargo test -p absrbrs lines::list`. Expected: compile errors (the types are missing).

- [ ] **Step 2: Implement**

Write the implementation at the top of `absrbrs/src/lines/list.rs`. Requirements:
- **Storage.** `molecule: Vec<u16>`, `isotopologue: Vec<u16>`, `columns: [Option<Vec<f64>>; 10]` indexed by `LineField as usize`. The seven required fields are always `Some`.
- **Validation in `from_columns`.** On any failure return `InvalidInput`, with a message naming the column and the first offending index. Check:
  - all columns, required and present optional, have the length of `nu`;
  - required fields are finite;
  - `nu > 0`, and `s_ref`, `gamma_air` and `gamma_self` are ≥ 0;
  - optional fields are finite or NaN; ±∞ is rejected.
- **Sorting.** Stable sort by `nu` (`f64::total_cmp`), applying one permutation to every column.
- **`value(field, i)`.** For optional fields, a missing column or a NaN entry resolves to the default: n_self → `value(NAir, i)`, δ_self → `value(DeltaAir, i)`, δ′_air → 0.0. `has_value` is false exactly in those cases. Both panic on `i >= len()`, like slice indexing; document this.
- **`select`.** Keeps rows that match all given filters: molecule equal; isotopologue in the slice; nu in the inclusive range. The result is still sorted.
- **`with_column`.** Validates `values` like `from_columns` (length, then the per-field rule) and returns a new list. Replacing `Nu` re-sorts. Replacing an optional field makes its column present.
- **Docs.** Every public item gets a doc comment. The `LineList` doc says that indices refer to the sorted order.

`LineField` docs: one line each, with the unit. For example, `Nu`: "line centre ν₀, cm⁻¹"; `SRef`: "intensity at 296 K, cm⁻¹/(molecule cm⁻²)"; `DeltapAir`: "temperature coefficient of the air shift, cm⁻¹ atm⁻¹ K⁻¹ (optional; default 0)".

- [ ] **Step 3: Run, check and commit**

Run:
```bash
cargo test -p absrbrs lines
cargo clippy --workspace --all-targets -- -D warnings
cargo +1.99.0 clippy --workspace --all-targets -- -D warnings
cargo fmt --all -- --check
RUSTDOCFLAGS="-D warnings" cargo doc -p absrbrs --no-deps --document-private-items
```
Expected: `9 passed` for the filter `lines` (3 from Task 1 and 6 new); everything else clean.

Commit:
```bash
git add absrbrs/src/lines
git commit -m "feat(lines): LineList (sorted SoA columns, optional fields with defaults, select, with_column)

Co-Authored-By: <your model> <noreply@anthropic.com>"
```

---

### Task 3: Readers (HITRAN `.par`, HAPI `.data`/`.header`)

**Files:**
- Modify: `Cargo.toml` (`serde_json` in `[workspace.dependencies]`), `absrbrs/Cargo.toml` (move `serde_json` to `[dependencies]`, keep the `float_roundtrip` feature)
- Create: `absrbrs/src/io/mod.rs`, `absrbrs/src/io/hitran.rs`, `absrbrs/src/io/hapi.rs`
- Modify: `absrbrs/src/lib.rs` (`pub mod io;`)
- Test: `absrbrs/tests/io_readers.rs`, `absrbrs/tests/local_hitran.rs`

**Interfaces:**
- Consumes: `LineList::from_columns`, `LineColumns`, and `Error::{Parse, Io}`.
- Produces:
  - `pub fn absrbrs::io::read_hitran_par(path: impl AsRef<Path>) -> Result<LineList>`;
  - `pub fn absrbrs::io::read_hapi(directory: impl AsRef<Path>, table: &str) -> Result<LineList>`, which reads `<directory>/<table>.header` and `<directory>/<table>.data`.

The record layout is the HITRAN 160-character format. Columns are zero-based byte offsets with widths, as in HAPI's `HITRAN_DEFAULT_HEADER`:

| Field | Bytes | Field | Bytes |
|---|---|---|---|
| molecule | 0..2 | gamma_self | 40..45 |
| isotopologue (1 char) | 2..3 | e_lower | 45..55 |
| nu | 3..15 | n_air | 55..59 |
| s_ref | 15..25 | delta_air | 59..67 |
| gamma_air | 35..40 | | |

The HAPI `.header` is JSON:
- `order` lists the fixed-width fields;
- `format` gives each field's `%<width>…` format;
- `position`, when present, gives the start offset of each field;
- `extra` lists the extra fields, in order;
- `extra_separator` is usually `,`.

A `.data` line is the fixed-width part, then the separator, then the extra values in `extra` order. `#` (or any unparsable token) means missing. The reader takes `nu`, `sw` (S_ref), `elower`, `gamma_air`, `gamma_self`, `n_air`, `delta_air`, `molec_id` and `local_iso_id` from the fixed part. When they appear in `extra`, it also takes `n_self`, `delta_self` and `deltap_air` (missing becomes NaN, meaning absent per line).

- [ ] **Step 1: Audit and add serde_json**

Run the usask-foundation:dependency-audit skill for `serde_json` 1 (MIT OR Apache-2.0; already a dev-dependency). Then:
- in the root `Cargo.toml`, add `serde_json = { version = "1", features = ["float_roundtrip"] }` to `[workspace.dependencies]`;
- in `absrbrs/Cargo.toml`, add `serde_json.workspace = true` to `[dependencies]` and remove `serde_json` from `[dev-dependencies]`. `serde` stays a dev-dependency.

- [ ] **Step 2: Write the failing tests**

Create `absrbrs/tests/io_readers.rs`. It writes synthetic records only; no HITRAN data is committed.
```rust
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
fn record(mol: u16, iso: char, nu: f64, s: f64, gair: f64, gself: f64, el: f64, nair: f64, dair: f64) -> String {
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
        record(2, '1', 6300.5, 1.234e-23, 0.0712, 0.093, 123.4567, 0.75, -0.0045),
        record(2, '0', 6200.25, 2.5e-25, 0.07, 0.09, 50.0, 0.7, -0.004),
        record(2, 'A', 6250.0, 3.0e-26, 0.07, 0.09, 50.0, 0.7, 0.0),
        record(2, '2', 6100.0, 0.0, 0.07, 0.09, 50.0, 0.7, 0.0), // S = 0: skipped
    ]
    .join("\n");
    std::fs::write(&path, text + "\n").unwrap();
    let lines = read_hitran_par(&path).unwrap();
    assert_eq!(lines.len(), 3);
    assert_eq!(lines.column(LineField::Nu).unwrap(), &[6200.25, 6250.0, 6300.5]);
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
    assert!(matches!(err, absrbrs::Error::Parse { line_number: 2, .. }), "{err}");
    std::fs::write(&path, format!("{good}\n{}\n", &good[..100])).unwrap();
    assert!(matches!(read_hitran_par(&path), Err(absrbrs::Error::Parse { line_number: 2, .. })));
}

#[test]
fn missing_file_is_an_io_error() {
    let err = read_hitran_par("/nonexistent/absrbrs/CO2.par").unwrap_err();
    assert!(matches!(err, absrbrs::Error::Io { kind: std::io::ErrorKind::NotFound, .. }));
}

/// A HAPI header with the default fixed part and extra columns, as `hapi.fetch` writes it.
fn hapi_header(extra: &[&str]) -> String {
    let order = ["molec_id", "local_iso_id", "nu", "sw", "a", "gamma_air", "gamma_self", "elower", "n_air", "delta_air"];
    let format = [("molec_id", "%2d"), ("local_iso_id", "%1d"), ("nu", "%12.6f"), ("sw", "%10.3E"), ("a", "%10.3E"),
        ("gamma_air", "%5.4f"), ("gamma_self", "%5.3f"), ("elower", "%10.4f"), ("n_air", "%4.2f"), ("delta_air", "%8.6f")];
    let position = [("molec_id", 0), ("local_iso_id", 2), ("nu", 3), ("sw", 15), ("a", 25), ("gamma_air", 35),
        ("gamma_self", 40), ("elower", 45), ("n_air", 55), ("delta_air", 59)];
    let extra_format: Vec<String> = extra.iter().map(|e| format!("\"{e}\": \"%9.6f\"")).collect();
    format!(
        "{{\"table_name\": \"T\", \"order\": {order:?}, \"format\": {{{}}}, \"position\": {{{}}}, \"extra\": {extra:?}, \"extra_format\": {{{}}}, \"extra_separator\": \",\"}}",
        format.iter().map(|(k, v)| format!("\"{k}\": \"{v}\"")).collect::<Vec<_>>().join(", "),
        position.iter().map(|(k, v)| format!("\"{k}\": {v}")).collect::<Vec<_>>().join(", "),
        extra_format.join(", ")
    )
}

#[test]
fn reads_hapi_caches_with_optional_extras() {
    let dir = tmp("hapi");
    std::fs::write(dir.join("CO.header"), hapi_header(&["gamma_sdv_0_air_296", "n_self", "delta_self"])).unwrap();
    let a = record(5, '1', 2140.0, 1e-19, 0.05, 0.06, 10.0, 0.7, -0.002);
    let b = record(5, '2', 2100.0, 1e-21, 0.05, 0.06, 10.0, 0.7, -0.002);
    std::fs::write(dir.join("CO.data"), format!("{a},0.0500,0.65,-0.001\n{b},#,#,0.003\n")).unwrap();
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
    std::fs::write(dir.join("X.data"), record(7, '1', 13100.0, 1e-23, 0.04, 0.04, 1.0, 0.7, 0.0) + "\n").unwrap();
    assert_eq!(read_hapi(&dir, "X").unwrap().len(), 1);
    std::fs::write(dir.join("Y.header"), "{not json").unwrap();
    std::fs::write(dir.join("Y.data"), "").unwrap();
    assert!(matches!(read_hapi(&dir, "Y"), Err(absrbrs::Error::Parse { .. })));
    assert!(matches!(read_hapi(&dir, "missing"), Err(absrbrs::Error::Io { .. })));
}
```
Also add a proptest in the same file: random strings of up to 300 bytes, mixing printable ASCII with the occasional non-ASCII byte, written as `.par` files, never panic. Allow only `Ok` or `Err(Parse | InvalidInput)`.
```rust
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
```

Create `absrbrs/tests/local_hitran.rs`, an opt-in test that is skipped unless the variable is set:
```rust
//! Opt-in: read the user's own HAPI caches (e.g. sasktran2's database/hitran_lines). Set
//! ABSRBRS_HITRAN_LINES_DIR to that directory; skipped otherwise (and never run in CI, since
//! HITRAN data is not redistributed).

use absrbrs::io::read_hapi;
use absrbrs::lines::{LineField, isotopologue};

#[test]
fn reads_local_hapi_caches() {
    let Ok(dir) = std::env::var("ABSRBRS_HITRAN_LINES_DIR") else {
        eprintln!("skipped: set ABSRBRS_HITRAN_LINES_DIR to a directory of HAPI caches");
        return;
    };
    let mut tables: Vec<_> = std::fs::read_dir(&dir).unwrap()
        .filter_map(|e| e.ok()?.path().file_name()?.to_str()?.strip_suffix(".header").map(String::from))
        .collect();
    tables.sort();
    for table in tables {
        let t = std::time::Instant::now();
        let lines = read_hapi(&dir, &table).unwrap();
        let nu = lines.column(LineField::Nu).unwrap();
        assert!(nu.windows(2).all(|w| w[0] <= w[1]));
        let unknown = (0..lines.len())
            .filter(|&i| isotopologue(lines.molecule()[i], lines.isotopologue()[i]).is_none())
            .count();
        eprintln!("{table}: {} lines in {:.2?}; {unknown} without isotopologue data", lines.len(), t.elapsed());
    }
}
```

Run: `cargo test -p absrbrs --test io_readers`. Expected: compile error, ``unresolved import `absrbrs::io` ``.

- [ ] **Step 3: Implement the readers**

Create `absrbrs/src/io/mod.rs`:
```rust
//! Line-file readers: HITRAN `.par` records and HAPI `.data`/`.header` caches (spec §5.6).

mod hapi;
mod hitran;

pub use hapi::read_hapi;
pub use hitran::read_hitran_par;
```
Add `pub mod io;` to `lib.rs`.

`absrbrs/src/io/hitran.rs` requirements:
- **Reading.** Open the file with `BufReader` and read line by line with `read_line` into one reused `String`. Map I/O failures to `Error::Io { path, kind, message }`, where `path` is the path as given (`path.display().to_string()`).
- **Per line (1-based `line_number`):**
  - strip a trailing `\r\n` or `\n`;
  - skip empty lines;
  - the record must be ASCII and at least 160 bytes long; otherwise `Parse { path, line_number, message: "record is N bytes; HITRAN records are 160" }`, or "contains non-ASCII bytes".
- **Fields.** Parse with `str::trim` and `parse::<f64>` (`parse::<u16>` for the molecule) from the byte ranges in the table above. A parse failure is `Parse`, with the message naming the field and the text, e.g. ``nu: cannot parse "abc.defgh"``.
- **Isotopologue character:** `'1'..='9'` gives 1–9, `'0'` gives 10, and `'A'..='Z'` gives 11–36. Any other character is `Parse`.
- **Skipping.** Records with `s_ref == 0.0` are skipped, as in sasktran2.
- **Result.** Collect the columns, then call `LineList::from_columns`. Its `InvalidInput`, for example a negative width, is returned as is.
- **Shared parser.** Expose the fixed-part parser as `pub(super) fn parse_fixed(record: &str, fields: &FixedLayout) -> std::result::Result<FixedRecord, String>`. The HAPI reader reuses it with its own offsets. `FixedLayout` gives the byte ranges of the nine fields; `FixedLayout::HITRAN` holds the defaults above.

`absrbrs/src/io/hapi.rs` requirements:
- **Header.** Read `<table>.header` and parse it with `serde_json::from_str::<serde_json::Value>`. Invalid JSON, or a missing `order`, is `Parse { path: <header path>, line_number: 1, .. }`.
- **Fixed layout.** Build a `FixedLayout` from `order`, `format` and `position`:
  - a field's width is the integer after `%` in its format;
  - its start is `position[field]` when `position` exists, and otherwise the running sum of the widths in `order`;
  - all nine fields must be present, or it is `Parse`, naming the field.
- **Extras.** `extra` (optional, default empty) and `extra_separator` (default `,`). For each of `n_self`, `delta_self` and `deltap_air` present in `extra`, note its index k.
- **Data lines.**
  - Strip the line ending and skip empty lines.
  - Parse the fixed part with the layout, using the same ASCII and length rules, with the length taken as the largest start + width.
  - If there are extras, split the line by the separator. Chunk 0 is the fixed part, and extra k is chunk k + 1. A missing chunk, or a token that does not parse as f64 (`#`), gives NaN.
  - Skip `s_ref == 0` as in the `.par` reader.
- **Result.** Columns for `n_self`, `delta_self` and `deltap_air` are `Some` exactly when the field is in `extra`.

- [ ] **Step 4: Run the tests and checks**

Run:
```bash
cargo test -p absrbrs --test io_readers
ABSRBRS_HITRAN_LINES_DIR="$HOME/Library/Application Support/sasktran2/database/hitran_lines" cargo test -p absrbrs --release --test local_hitran -- --nocapture
cargo test -p absrbrs
cargo clippy --workspace --all-targets -- -D warnings && cargo +1.99.0 clippy --workspace --all-targets -- -D warnings
cargo fmt --all -- --check
```
Expected:
- `io_readers`: `6 passed`, including the proptest.
- `local_hitran`: one line per table, e.g. `CO: 5381 lines in …; 0 without isotopologue data`. The large tables (CH₄, CO₂, H₂O, 100–170 MB) should read in a few seconds each in release mode. Report the timings.
- Everything else passes and is clean.

If the local test isn't available on the machine running it, it prints `skipped` and passes.

- [ ] **Step 5: Commit**

```bash
git add Cargo.toml Cargo.lock absrbrs/Cargo.toml absrbrs/src/lib.rs absrbrs/src/io absrbrs/tests/io_readers.rs absrbrs/tests/local_hitran.rs
git commit -m "feat(io): HITRAN .par and HAPI .data/.header readers

Dependency audit: serde_json 1 (MIT OR Apache-2.0) ...

Co-Authored-By: <your model> <noreply@anthropic.com>"
```

---

### Task 4: TIPS-2021 partition functions

**Files:**
- Create: `tools/tips2021.py`, `absrbrs/data/tips2021.bin` and `absrbrs/data/tips2021.json` (generated), `absrbrs/LICENSES/TIPS-2021-MIT.txt`
- Create: `absrbrs/src/partition/{mod.rs, hermite.rs, tips2021.rs, tabulated.rs}`
- Modify: `absrbrs/src/lib.rs` (`pub mod partition;`), `pyproject.toml` (`license-files`; `[tool.uv] cache-keys` gains `{ file = "absrbrs/data/**" }` and `{ file = "rust-toolchain.toml" }`), `THIRD_PARTY.md` (TIPS row: files)
- Test: unit tests in the partition modules, `tests/test_data_files.py` (Python)

**Interfaces:**
- Produces, in `absrbrs::partition`:
  - `pub struct PartitionValue { pub q: f64, pub dlnq_dt: f64 }`, deriving `Clone, Copy, Debug, PartialEq`;
  - `pub trait PartitionFunction: Send + Sync { fn evaluate(&self, molecule: u16, isotopologue: u16, temperature: f64) -> Result<PartitionValue>; }`;
  - `pub struct Tips2021`, with `pub fn tips2021() -> &'static Tips2021`, `Tips2021::temperature_range(&self, molecule: u16, isotopologue: u16) -> Option<(f64, f64)>` and `Tips2021::isotopologues(&self) -> Vec<(u16, u16)>`;
  - `pub struct TabulatedPartition`, deriving `Clone, Debug, Default`, with `new()` and `insert(&mut self, molecule: u16, isotopologue: u16, temperature: &[f64], q: &[f64]) -> Result<()>`;
  - both implement `PartitionFunction`.
- Errors:
  - `MissingPartitionFunction { molecule, isotopologue }` when there is no table;
  - `OutsideTable { molecule, isotopologue, temperature }` when T is outside [first node, last node], or is NaN.
- Embedded file format (all little-endian):
  - the 16-byte magic `b"ABSRBRS-TIPS2021"`;
  - `u32` count;
  - count entries of 16 bytes: `u16` molecule, `u16` isotopologue, `u32` n_points, `f32` t_first, `f32` t_step;
  - then the f32 Q values for each entry in turn.

- [ ] **Step 1: Write the converter and generate the data**

Create `tools/tips2021.py`:
```python
# /// script
# requires-python = ">=3.11"
# dependencies = []
# ///
"""Build absrbrs/data/tips2021.bin from the TIPS-2021 tables (Gamache et al. 2021).

Source: Zenodo record 10.5281/zenodo.4708099 (MIT licence, Copyright (C) 2021 Robert Gamache),
one `<molecule>_<isotopologue>.QTpy` file per isotopologue: a pickled dict of strings, temperature
(K, 1 K steps from 1 K) to Q. Each file's MD5 is checked against Zenodo's record, its opcodes are
inspected (no globals or reductions allowed), and it is loaded with an unpickler that refuses
all globals. Q is stored as f32 (≤ 6e-8 relative rounding, at or below the tables' 6-9
significant digits). Run with `uv run tools/tips2021.py`; the outputs are committed.
"""

from __future__ import annotations

import hashlib
import io
import json
import pathlib
import pickle
import pickletools
import struct
import urllib.request

RECORD = "4708099"
API = f"https://zenodo.org/api/records/{RECORD}"
UA = {"User-Agent": "absrbrs-tools"}
DATA = pathlib.Path(__file__).resolve().parents[1] / "absrbrs" / "data"
MAGIC = b"ABSRBRS-TIPS2021"
FORBIDDEN = {"GLOBAL", "STACK_GLOBAL", "REDUCE", "INST", "OBJ", "NEWOBJ", "NEWOBJ_EX", "BUILD", "EXT1", "EXT2", "EXT4"}


class Refuse(pickle.Unpickler):
    def find_class(self, module, name):  # noqa: ARG002
        raise pickle.UnpicklingError(f"global {module}.{name} refused")


def fetch(url: str) -> bytes:
    with urllib.request.urlopen(urllib.request.Request(url, headers=UA)) as r:
        return r.read()


def load_table(data: bytes) -> list[tuple[float, float]]:
    ops = {op.name for op, _, _ in pickletools.genops(data)}
    if ops & FORBIDDEN:
        raise ValueError(f"unexpected pickle opcodes {sorted(ops & FORBIDDEN)}")
    table = Refuse(io.BytesIO(data)).load()
    if not isinstance(table, dict) or not all(isinstance(k, str) and isinstance(v, str) for k, v in table.items()):
        raise ValueError("expected a dict of str to str")
    return sorted((float(k), float(v)) for k, v in table.items())


def main() -> None:
    record = json.loads(fetch(API))
    files = {
        f["key"]: f
        for f in record["files"]
        if f["key"].endswith(".QTpy") and f["key"] != "QT_HITRAN.QTpy"
    }
    entries, values, sources = [], [], {}
    for key in sorted(files, key=lambda k: tuple(int(p) for p in k[:-5].split("_"))):
        meta = files[key]
        data = fetch(meta["links"]["self"])
        md5 = hashlib.md5(data).hexdigest()  # Zenodo publishes MD5
        if f"md5:{md5}" != meta["checksum"]:
            raise ValueError(f"{key}: MD5 {md5} does not match Zenodo's {meta['checksum']}")
        rows = load_table(data)
        temps = [t for t, _ in rows]
        if temps != [1.0 + i for i in range(len(temps))]:
            raise ValueError(f"{key}: temperatures are not 1, 2, …, N K")
        mol, iso = (int(p) for p in key[:-5].split("_"))
        entries.append((mol, iso, len(rows), 1.0, 1.0))
        values.extend(q for _, q in rows)
        sources[key] = {"md5": md5, "sha256": hashlib.sha256(data).hexdigest(), "points": len(rows)}
    blob = bytearray(MAGIC)
    blob += struct.pack("<I", len(entries))
    for mol, iso, n, t0, dt in entries:
        blob += struct.pack("<HHIff", mol, iso, n, t0, dt)
    blob += struct.pack(f"<{len(values)}f", *values)
    DATA.mkdir(parents=True, exist_ok=True)
    (DATA / "tips2021.bin").write_bytes(bytes(blob))
    manifest = {
        "source": f"https://doi.org/10.5281/zenodo.{RECORD}",
        "licence": "MIT, Copyright (C) 2021 Robert Gamache (see absrbrs/LICENSES/TIPS-2021-MIT.txt)",
        "isotopologues": len(entries),
        "points": len(values),
        "sha256": hashlib.sha256(bytes(blob)).hexdigest(),
        "files": sources,
    }
    (DATA / "tips2021.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"{len(entries)} isotopologues, {len(values)} points, {len(blob)} bytes")


if __name__ == "__main__":
    main()
```
Run: `uv run tools/tips2021.py`. It downloads about 36 MB from Zenodo.

Expected: `181 isotopologues, <N> points, <bytes> bytes`, with about 1 million points and about 4.2 MB. If an MD5 does not match, stop and report. Do not change the check.

Create `absrbrs/LICENSES/TIPS-2021-MIT.txt`: the heading `TIPS-2021 total internal partition sums (absrbrs/data/tips2021.bin is converted from the Zenodo record 10.5281/zenodo.4708099)`, a blank line, then the copyright and MIT text from the header of `TIPS_2021.py` in the same Zenodo record. Fetch that file with the record's file link and copy the notice exactly.

Add the file to `pyproject.toml` `license-files`. In `THIRD_PARTY.md`'s TIPS-2021 row, add: "Embedded as `absrbrs/data/tips2021.bin` (f32), converted by `tools/tips2021.py` from the `.QTpy` files (MD5 and SHA-256 per file in `absrbrs/data/tips2021.json`); notice `absrbrs/LICENSES/TIPS-2021-MIT.txt`."

- [ ] **Step 2: Write the failing tests**

Create `tests/test_data_files.py` (root pytest):
```python
"""The embedded TIPS-2021 file matches its manifest."""

import hashlib
import json
import struct
from pathlib import Path

DATA = Path(__file__).resolve().parents[1] / "absrbrs" / "data"


def test_tips2021_blob_matches_its_manifest():
    blob = (DATA / "tips2021.bin").read_bytes()
    manifest = json.loads((DATA / "tips2021.json").read_text("utf-8"))
    assert hashlib.sha256(blob).hexdigest() == manifest["sha256"]
    assert blob[:16] == b"ABSRBRS-TIPS2021"
    (count,) = struct.unpack_from("<I", blob, 16)
    assert count == manifest["isotopologues"] == 181
    points = sum(struct.unpack_from("<HHIff", blob, 20 + 16 * k)[2] for k in range(count))
    assert points == manifest["points"]
    assert len(blob) == 20 + 16 * count + 4 * points
```

Create `absrbrs/src/partition/mod.rs`, with the tests at the end:
```rust
//! Partition functions Q(T) and d ln Q/dT (spec §5.5).

mod hermite;
mod tabulated;
mod tips2021;

pub use tabulated::TabulatedPartition;
pub use tips2021::{Tips2021, tips2021};

use crate::Result;

/// Q(T) and its logarithmic derivative d ln Q/dT (K⁻¹), consistent with each other.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PartitionValue {
    /// The total internal partition sum Q(T).
    pub q: f64,
    /// d ln Q/dT, K⁻¹.
    pub dlnq_dt: f64,
}

/// A source of partition functions, by HITRAN molecule and isotopologue number.
pub trait PartitionFunction: Send + Sync {
    /// Q and d ln Q/dT at `temperature` (K). Errors: `MissingPartitionFunction` if there is no
    /// table for the isotopologue, `OutsideTable` if the temperature is outside it (or NaN).
    fn evaluate(&self, molecule: u16, isotopologue: u16, temperature: f64) -> Result<PartitionValue>;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Error;

    #[test]
    fn tips_reproduces_the_source_at_nodes() {
        // Values from the Zenodo .QTpy files (1_1: Q(296) = 174.58130; 2_1: 286.09382;
        // 7_1: 215.734504; 3_1: 3474.9978), within f32 rounding.
        for (mol, iso, q296) in [(1, 1, 174.58130), (2, 1, 286.09382), (7, 1, 215.734504), (3, 1, 3474.9978)] {
            let v = tips2021().evaluate(mol, iso, 296.0).unwrap();
            assert!((v.q - q296).abs() <= 6e-8 * q296, "({mol}, {iso}): {} vs {q296}", v.q);
        }
        assert_eq!(tips2021().evaluate(1, 1, 1.0).unwrap().q, 1.0);
    }

    #[test]
    fn tips_derivative_is_that_of_the_interpolant() {
        for (mol, iso) in [(1, 1), (2, 1), (7, 1), (3, 1), (2, 10)] {
            for t in [3.7, 19.5, 150.25, 296.0, 299.9, 777.7] {
                let h = 1e-4;
                let lp = tips2021().evaluate(mol, iso, t + h).unwrap().q.ln();
                let lm = tips2021().evaluate(mol, iso, t - h).unwrap().q.ln();
                let d = tips2021().evaluate(mol, iso, t).unwrap().dlnq_dt;
                assert!(((lp - lm) / (2.0 * h) - d).abs() <= 1e-6 * d.abs().max(1e-3), "({mol}, {iso}) at {t}");
            }
        }
    }

    #[test]
    fn tips_coverage_and_errors() {
        let isos = tips2021().isotopologues();
        assert_eq!(isos.len(), 181);
        assert_eq!(tips2021().temperature_range(1, 1), Some((1.0, 5000.0)));
        assert_eq!(tips2021().temperature_range(3, 1), Some((1.0, 1000.0)));
        assert!(matches!(tips2021().evaluate(3, 1, 1000.5), Err(Error::OutsideTable { molecule: 3, isotopologue: 1, .. })));
        assert!(matches!(tips2021().evaluate(1, 1, 0.5), Err(Error::OutsideTable { .. })));
        assert!(matches!(tips2021().evaluate(1, 1, f64::NAN), Err(Error::OutsideTable { .. })));
        assert!(matches!(tips2021().evaluate(10, 3, 296.0), Err(Error::MissingPartitionFunction { molecule: 10, isotopologue: 3 })));
    }

    #[test]
    fn tabulated_tables() {
        let mut p = TabulatedPartition::new();
        // Q = T^1.5 exactly at the nodes; nonuniform grid.
        let t: Vec<f64> = vec![100.0, 150.0, 220.0, 300.0, 400.0];
        let q: Vec<f64> = t.iter().map(|t: &f64| t.powf(1.5)).collect();
        p.insert(2, 1, &t, &q).unwrap();
        let v = p.evaluate(2, 1, 220.0).unwrap();
        assert!((v.q - 220f64.powf(1.5)).abs() <= 1e-12 * v.q);
        let v = p.evaluate(2, 1, 250.0).unwrap();
        assert!((v.q / 250f64.powf(1.5) - 1.0).abs() < 1e-3);
        assert!((v.dlnq_dt - 1.5 / 250.0).abs() < 1e-4);
        assert!(matches!(p.evaluate(2, 2, 250.0), Err(Error::MissingPartitionFunction { .. })));
        assert!(matches!(p.evaluate(2, 1, 401.0), Err(Error::OutsideTable { .. })));
        assert!(p.insert(2, 3, &[1.0, 1.0], &[1.0, 2.0]).is_err()); // not increasing
        assert!(p.insert(2, 3, &[1.0, 2.0], &[1.0, -2.0]).is_err()); // Q must be > 0
        assert!(p.insert(2, 3, &[1.0], &[1.0]).is_err()); // at least two points
        assert!(p.insert(2, 3, &[1.0, 2.0], &[1.0]).is_err()); // lengths
    }
}
```
Add `pub mod partition;` to `lib.rs`.

Run: `cargo test -p absrbrs partition` and `uv run --locked pytest -q tests/test_data_files.py`. Expected: Rust compile errors, because the modules are missing; the Python test passes.

- [ ] **Step 3: Implement**

Create `absrbrs/src/partition/hermite.rs`. It holds the shared interpolant: cubic Hermite in ln Q with node slopes from numpy's second-order `gradient` formula, using one-sided first-order slopes at the two ends.
```rust
//! Cubic Hermite interpolation of ln Q(T), with node slopes from the second-order centred
//! difference (numpy's `gradient`; first order at the ends). The derivative returned is that of
//! the interpolant, so Q and d ln Q/dT are exactly consistent.

/// A table of ln Q at increasing temperatures.
#[derive(Clone, Debug)]
pub(crate) struct Table {
    pub(crate) t: Vec<f64>,
    pub(crate) ln_q: Vec<f64>,
    slope: Vec<f64>,
}

impl Table {
    /// `t` strictly increasing with at least two points, `ln_q` the same length (checked by callers).
    pub(crate) fn new(t: Vec<f64>, ln_q: Vec<f64>) -> Table {
        let n = t.len();
        let mut slope = vec![0.0; n];
        slope[0] = (ln_q[1] - ln_q[0]) / (t[1] - t[0]);
        slope[n - 1] = (ln_q[n - 1] - ln_q[n - 2]) / (t[n - 1] - t[n - 2]);
        for i in 1..n - 1 {
            let (hl, hr) = (t[i] - t[i - 1], t[i + 1] - t[i]);
            slope[i] = (hl * hl * ln_q[i + 1] - hr * hr * ln_q[i - 1] + (hr * hr - hl * hl) * ln_q[i])
                / (hl * hr * (hl + hr));
        }
        Table { t, ln_q, slope }
    }

    /// The first and last temperatures.
    pub(crate) fn range(&self) -> (f64, f64) {
        (self.t[0], self.t[self.t.len() - 1])
    }

    /// (ln Q, d ln Q/dT) at `temperature`, which must lie in `range()` (checked by callers).
    pub(crate) fn evaluate(&self, temperature: f64) -> (f64, f64) {
        let n = self.t.len();
        let i = self.t.partition_point(|&t| t <= temperature).clamp(1, n - 1) - 1;
        let h = self.t[i + 1] - self.t[i];
        let s = (temperature - self.t[i]) / h;
        let (y0, y1, m0, m1) = (self.ln_q[i], self.ln_q[i + 1], self.slope[i] * h, self.slope[i + 1] * h);
        let s2 = s * s;
        let s3 = s2 * s;
        let value = (2.0 * s3 - 3.0 * s2 + 1.0) * y0 + (s3 - 2.0 * s2 + s) * m0 + (-2.0 * s3 + 3.0 * s2) * y1 + (s3 - s2) * m1;
        let deriv = ((6.0 * s2 - 6.0 * s) * y0 + (3.0 * s2 - 4.0 * s + 1.0) * m0 + (-6.0 * s2 + 6.0 * s) * y1 + (3.0 * s2 - 2.0 * s) * m1) / h;
        (value, deriv)
    }
}
```
For uniform grids the slope formula reduces to the centred difference. Q at a node is `exp(ln Q)`, which equals the stored value to about 1 ulp.

Create `absrbrs/src/partition/tips2021.rs`:
- `static BLOB: &[u8] = include_bytes!("../../data/tips2021.bin");`
- `pub struct Tips2021 { index: BTreeMap<(u16, u16), Entry> }`, where `Entry { offset: usize, n: usize, t_first: f64, t_step: f64, table: OnceLock<Table> }`. Each table is decoded on first use: ln of every f32 value as f64, with t_k = t_first + k·t_step.
- `pub fn tips2021() -> &'static Tips2021` uses a `static TIPS: OnceLock<Tips2021>` and parses the header once. The blob is compile-time data that the Python test checks, so parse failures may use `expect` with a message saying the embedded file is malformed. Document this in a comment.
- `evaluate`:
  - with no entry, return `MissingPartitionFunction`;
  - if `!(t_first <= T && T <= t_last)`, return `OutsideTable`;
  - otherwise `(ln_q, d) = table.evaluate(T)` and return `PartitionValue { q: ln_q.exp(), dlnq_dt: d }`.
- `temperature_range` and `isotopologues()` (sorted) come from the index.
- Doc comment on `Tips2021`: "Total internal partition sums from TIPS-2021 (Gamache et al. 2021; MIT licence, Copyright (C) 2021 Robert Gamache), for 181 HITRAN isotopologues at 1 K steps from 1 K to each table's maximum (1000–7500 K), stored as f32 (≤ 6e-8 relative)."

Create `absrbrs/src/partition/tabulated.rs`: `TabulatedPartition { tables: BTreeMap<(u16, u16), Table> }`. `insert` validates:
- the lengths are equal and at least 2;
- `temperature` is finite, strictly increasing and > 0;
- `q` is finite and > 0.

It returns `InvalidInput` with a clear message on failure, and otherwise stores `Table::new(t, q.ln())`, replacing any previous table. `evaluate` behaves as in TIPS.

In `pyproject.toml` `[tool.uv] cache-keys`, add `{ file = "absrbrs/data/**" }` and `{ file = "rust-toolchain.toml" }`. These were a phase 0 review note, and the embedded data must trigger a rebuild.

- [ ] **Step 4: Run, check package size, commit**

Run:
```bash
cargo test -p absrbrs partition
uv run --locked pytest -q tests/test_data_files.py
cargo clippy --workspace --all-targets -- -D warnings && cargo +1.99.0 clippy --workspace --all-targets -- -D warnings
cargo fmt --all -- --check
cargo package -p absrbrs --allow-dirty --list | grep -E "data/tips2021|LICENSES"
cargo package -p absrbrs --allow-dirty && ls -l target/package/absrbrs-*.crate
```
Expected:
- 4 Rust partition tests pass;
- the Python test passes;
- the package lists `data/tips2021.bin`, `data/tips2021.json` and the three `LICENSES` files;
- the `.crate` is well under 10 MB. Report its size.

Commit:
```bash
git add tools/tips2021.py absrbrs/data absrbrs/LICENSES/TIPS-2021-MIT.txt absrbrs/src/lib.rs absrbrs/src/partition tests/test_data_files.py pyproject.toml THIRD_PARTY.md
git commit -m "feat(partition): TIPS-2021 (181 isotopologues, 1 K, f32) and user tables, cubic Hermite in ln Q

Co-Authored-By: <your model> <noreply@anthropic.com>"
```

---

### Task 5: Line-state layer

**Files:**
- Create: `absrbrs/src/state/{mod.rs, real.rs, jacobian.rs, dual.rs}`
- Modify: `absrbrs/src/lib.rs` (`mod state;`, crate-private)
- Create: `reference/absrbrs_reference/generators/line_state.py`; modify `generators/__init__.py`; generate `reference/data/line_state.jsonl`
- Test: unit tests in `state/`: the dual-number checks, and the Arb comparison (a unit test, since integration tests cannot reach crate-private items; it reads the reference file through `crate::test_data::reference_file`)

**Interfaces:**
- Consumes: `crate::lines::{LineList, LineField}`.
- Produces (crate-private, in `crate::state`):
  - `pub(crate) trait Real`;
  - `pub(crate) struct Condition<R> { pub temperature: R, pub pressure: R, pub vmr: R }`;
  - `pub(crate) struct LineParameters<R> { pub nu, s_ref, e_lower, gamma_air, gamma_self, n_air: R, pub n_self: Option<R>, pub delta_air: R, pub deltap_air: R, pub delta_self: Option<R> }`, plus `LineParameters::from_list(&LineList, i) -> LineParameters<f64>`, where n_self and δ_self are `None` when absent;
  - `pub(crate) struct LineShape<R> { pub amplitude: R, pub centre: R, pub doppler_width: R, pub y: R }`;
  - `pub(crate) fn line_shape<R: Real>(c: &Condition<R>, l: &LineParameters<R>, log_q_ratio: R, mass: f64) -> LineShape<R>`, where `log_q_ratio` = ln(Q(296)/Q(T)) and `mass` is in g/mol;
  - `pub(crate) struct PartitionRatio { pub log_q_ratio: f64, pub dlnq_dt: f64 }`;
  - `pub(crate) struct StateJacobian { pub temperature, pressure, vmr: LineShape<f64> }`;
  - `pub(crate) fn line_shape_state_jacobian(c: &Condition<f64>, l: &LineParameters<f64>, p: PartitionRatio, mass: f64) -> (LineShape<f64>, StateJacobian)`;
  - `pub(crate) fn line_shape_line_jacobian(c: &Condition<f64>, l: &LineParameters<f64>, p: PartitionRatio, mass: f64) -> (LineShape<f64>, [LineShape<f64>; 10])`, indexed by `LineField as usize`.

Formulas (spec §5.2). Pressures are in atm: p_atm = p/101325, p_self = VMR·p_atm, p_air = (1 − VMR)·p_atm.
- **Intensity factor:** f = S(T)/S_ref = exp(log_q_ratio) · exp(−c₂E″(1/T − 1/T_ref)) · expm1(−c₂ν₀/T)/expm1(−c₂ν₀/T_ref).
- **Doppler width:** σ_D = (ν₀/c) √(2k_B T / m), with m = mass·10⁻³/N_A kg.
- **Lorentz width:** γ_L = (T_ref/T)^{n_air} γ_air p_air + (T_ref/T)^{n_self} γ_self p_self.
- **Shifted centre:** ν₀′ = ν₀ + (δ_air + δ′_air (T − T_ref)) p_air + δ_self p_self.
- **Shape quantities:** A = 10⁻⁴ S_ref f/(√π σ_D) and y = γ_L/σ_D. The cross section is σ(ν) = A·K(x, y), in m²/molecule.

- [ ] **Step 1: Write the Real trait and the dual-number test type**

Create `absrbrs/src/state/real.rs`:
```rust
//! The scalar arithmetic the line-state formulas need: f64 in production, dual numbers in tests.

use std::ops::{Add, Div, Mul, Neg, Sub};

pub(crate) trait Real:
    Copy + Add<Output = Self> + Sub<Output = Self> + Mul<Output = Self> + Div<Output = Self> + Neg<Output = Self>
{
    fn from_f64(v: f64) -> Self;
    fn exp(self) -> Self;
    fn exp_m1(self) -> Self;
    fn ln(self) -> Self;
    fn sqrt(self) -> Self;
    fn powf(self, e: Self) -> Self;
}

impl Real for f64 {
    fn from_f64(v: f64) -> f64 { v }
    fn exp(self) -> f64 { f64::exp(self) }
    fn exp_m1(self) -> f64 { f64::exp_m1(self) }
    fn ln(self) -> f64 { f64::ln(self) }
    fn sqrt(self) -> f64 { f64::sqrt(self) }
    fn powf(self, e: f64) -> f64 { f64::powf(self, e) }
}
```
Create `absrbrs/src/state/dual.rs`, marked `#[cfg(test)]`: `Dual<const N: usize> { v: f64, d: [f64; N] }` with:
- `Dual::constant(v)` and `Dual::variable(v, k)`, which sets d[k] = 1;
- the operations `Add`, `Sub`, `Mul`, `Div` and `Neg`, using the usual rules;
- `Real`:
  - exp: (eᵛ, eᵛ d);
  - exp_m1: (expm1 v, eᵛ d);
  - ln: (ln v, d/v);
  - sqrt: (√v, d/(2√v));
  - powf(e): w = vᵉ, d_w = w (d_e ln v + e d/v).

Keep it to about 80 lines.

- [ ] **Step 2: Write the generic formula and its dual-number tests**

Create `absrbrs/src/state/mod.rs` with the structs from the interfaces block, the constants below, and `line_shape`:
```rust
//! The line-state layer (spec §4.2, §5): from a condition and a line to the shape quantities
//! q = (A, ν₀′, σ_D, y), with closed-form derivatives in `jacobian`.

mod jacobian;
mod real;
#[cfg(test)]
mod dual;

pub(crate) use jacobian::{line_shape_line_jacobian, line_shape_state_jacobian, PartitionRatio, StateJacobian};
pub(crate) use real::Real;

pub(crate) const T_REF: f64 = 296.0;
pub(crate) const ATM: f64 = 101_325.0;
/// Second radiation constant hc/k_B, cm K (exact SI constants; nearest double).
pub(crate) const C2: f64 = 1.438_776_877_503_933_8;
pub(crate) const C_LIGHT: f64 = 299_792_458.0;
pub(crate) const K_B: f64 = 1.380_649e-23;
pub(crate) const N_A: f64 = 6.022_140_76e23;
const FRAC_1E4_SQRT_PI: f64 = 1e-4 * crate::faddeeva::FRAC_1_SQRT_PI;

/// q for one line at one condition (see the module docs for the formulas).
pub(crate) fn line_shape<R: Real>(c: &Condition<R>, l: &LineParameters<R>, log_q_ratio: R, mass: f64) -> LineShape<R> {
    let f = R::from_f64;
    let p = c.pressure / f(ATM);
    let p_self = c.vmr * p;
    let p_air = (f(1.0) - c.vmr) * p;
    let t = c.temperature;
    let tr = f(T_REF);
    let stim = (-(f(C2) * l.nu / t)).exp_m1() / (-(f(C2) * l.nu / tr)).exp_m1();
    let boltz = (-(f(C2) * l.e_lower) * (f(1.0) / t - f(1.0) / tr)).exp();
    let intensity = log_q_ratio.exp() * boltz * stim;
    let sigma_d = l.nu / f(C_LIGHT) * (f(2.0 * K_B * N_A * 1e3 / mass) * t).sqrt();
    let ratio = tr / t;
    let n_self = l.n_self.unwrap_or(l.n_air);
    let delta_self = l.delta_self.unwrap_or(l.delta_air);
    let gamma = ratio.powf(l.n_air) * l.gamma_air * p_air + ratio.powf(n_self) * l.gamma_self * p_self;
    LineShape {
        amplitude: f(FRAC_1E4_SQRT_PI) * l.s_ref * intensity / sigma_d,
        centre: l.nu + (l.delta_air + l.deltap_air * (t - tr)) * p_air + delta_self * p_self,
        doppler_width: sigma_d,
        y: gamma / sigma_d,
    }
}
```
Here `2.0 * K_B * N_A * 1e3 / mass` is 2k/m with m = mass·10⁻³/N_A.

In `state/mod.rs`, add a test module `dual_checks`:
- Build 13-component duals: T, p, VMR, then the ten `LineField`s in `LineField::ALL` order.
- `log_q_ratio` is a dual with value L and d[0] = −dlnq_dt (because ln(Q(T_ref)/Q(T)) falls with ln Q(T)).
- For each case:
  - run `line_shape` on the duals;
  - run the closed forms (`line_shape_state_jacobian`, `line_shape_line_jacobian`) with `PartitionRatio { log_q_ratio: L, dlnq_dt }`;
  - compare every component: |closed − dual| ≤ 1e-13·(|closed| + |dual|) + 1e-300.
  - The values must also equal the f64 `line_shape` bitwise.
- Cases, every combination of:
  - lines: (ν₀, S_ref, E″, γ_air, γ_self, n_air, n_self, δ_air, δ′_air, δ_self) =
    - (13100.5, 7.4e-24, 1200.0, 0.045, 0.047, 0.72, Some 0.70, −0.008, Some 1e-5, Some −0.007)
    - (6300.25, 1.6e-23, 300.0, 0.07, 0.09, 0.75, None, −0.004, 0.0, None)
    - (0.75, 1e-24, 45.0, 0.1, 0.4, 0.68, None, 0.0, 0.0, None), a microwave line, where c₂ν₀/T is small
    - (2140.0, 4e-19, 8000.0, 0.05, 0.06, 0.7, Some 0.6, −0.002, 0.0, Some 0.003)
  - conditions: T ∈ {150, 296, 1500}; p ∈ {0, 1000, 101325}; VMR ∈ {0, 0.01, 1}.
- Tie semantics, for a line with n_self and δ_self absent:
  - every closed-form component except `NSelf` and `DeltaSelf` must match the dual run where they are `None`;
  - the `NSelf` and `DeltaSelf` components must match a dual run where `n_self = Some(n_air value, seeded as NSelf)` and `delta_self = Some(delta_air value, seeded as DeltaSelf)`. That is the partial taken at the default.

Run: `cargo test -p absrbrs state`. Expected: compile errors, because `jacobian` is missing.

- [ ] **Step 3: Write the closed forms**

Create `absrbrs/src/state/jacobian.rs`. Compute `q = line_shape::<f64>(...)` first; returning that guarantees bitwise-equal values. Then compute these intermediates:
- p_atm, p_air, p_self, ratio = T_ref/T, ln_ratio = ln(ratio);
- g_air = ratio^{n_air}·γ_air and g_self = ratio^{n_self,eff}·γ_self;
- shift_air = δ_air + δ′_air (T − T_ref);
- u = c₂ν₀/T and u_r = c₂ν₀/T_ref;
- `unit_amplitude` = 10⁻⁴ f/(√π σ_D), the amplitude per unit S_ref. Compute it directly rather than as A/S_ref, so that S_ref = 0 works.

For every θ, `∂y = (∂γ_L − y ∂σ_D)/σ_D`.

**State:**

| θ | ∂A | ∂σ_D | ∂γ_L | ∂ν₀′ |
|---|---|---|---|---|
| T | A·(−dlnq_dt + c₂E″/T² − (u/T)/expm1(u) − 1/(2T)) | σ_D/(2T) | −(n_air/T) g_air p_air − (n_self,eff/T) g_self p_self | δ′_air p_air |
| p (Pa) | 0 | 0 | (g_air (1 − VMR) + g_self VMR)/101325 | (shift_air (1 − VMR) + δ_self,eff VMR)/101325 |
| VMR | 0 | 0 | p_atm (g_self − g_air) | p_atm (δ_self,eff − shift_air) |

**Line fields.** Only nonzero entries are listed; the rest are 0. "Tied" means the optional field is absent (`None`).

| θ | Effect |
|---|---|
| Nu | ∂ν₀′ = 1; ∂σ_D = σ_D/ν₀; ∂A = A·((c₂/T)/expm1(u) − (c₂/T_ref)/expm1(u_r) − 1/ν₀) |
| SRef | ∂A = unit_amplitude |
| ELower | ∂A = −A·c₂(1/T − 1/T_ref) |
| GammaAir | ∂γ_L = ratio^{n_air} p_air |
| GammaSelf | ∂γ_L = ratio^{n_self,eff} p_self |
| NAir | ∂γ_L = ln_ratio·g_air p_air, plus ln_ratio·g_self p_self if n_self is tied |
| NSelf | ∂γ_L = ln_ratio·g_self p_self |
| DeltaAir | ∂ν₀′ = p_air, plus p_self if δ_self is tied |
| DeltapAir | ∂ν₀′ = (T − T_ref) p_air |
| DeltaSelf | ∂ν₀′ = p_self |

In the `Nu` and `T` rows, `(u/T)/expm1(u)` and `(c₂/T)/expm1(u)` must stay finite for large u, where expm1 → ∞ and the result → 0. Write them as `u / T / u.exp_m1()`, which gives 0 when the denominator is ∞, and never as `exp(-u)/(1 - exp(-u))`.

Run: `cargo test -p absrbrs state -- --nocapture`. Expected: every dual check passes.

- [ ] **Step 4: Arb reference values**

Create `reference/absrbrs_reference/generators/line_state.py`. It evaluates the formulas above in Arb at 256 bits, an implementation independent of the Rust code:
- c₂ is the exact rational hc/k_B·100 (cm K) from h = 6.62607015e-34, c = 299792458 and k = 1.380649e-23;
- m = mass·10⁻³/N_A, exactly;
- every input is a decimal literal parsed by Python (exact doubles, with no libm involved).

Cases: the same lines and conditions as the dual checks, plus T = 50 K and p = 1e6 Pa. Inputs per record:
- `temperature`, `pressure`, `vmr`;
- the line fields, with `n_self` and `delta_self` as `null` when absent;
- `log_q_ratio` (choose 0.0, 0.37 and −1.2, one per line, cycling) and `mass` (32.0, 44.0, 18.0, 28.0).

Outputs: `amplitude`, `centre`, `doppler_width`, `y`, using `records.record` with the label `certified`. Register it as `"line_state": line_state.generate` and run `uv run --locked --directory reference python -m absrbrs_reference generate --only line_state`.

In `state/mod.rs`, add a unit test `matches_arb_reference`. It reads `line_state.jsonl` through `crate::test_data::reference_file` and asserts relative error ≤ 1e-13 for each of the four outputs.

- [ ] **Step 5: Run, check, commit**

Run:
```bash
cargo test -p absrbrs
uv run --locked --directory reference python -m pytest -q
cargo clippy --workspace --all-targets -- -D warnings && cargo +1.99.0 clippy --workspace --all-targets -- -D warnings
cargo fmt --all -- --check
uv run --locked ruff check . && uv run --locked ruff format --check .
```
Expected: all pass. Report the largest relative error against Arb.

`state` has no consumer outside tests until phase 1c. Mark its items `#[allow(dead_code)]` with the comment `// used by the engine (phase 1c)`, rather than making them public.

Commit:
```bash
git add absrbrs/src/lib.rs absrbrs/src/state reference
git commit -m "feat(state): line-state layer with closed-form state and line-parameter Jacobians

Co-Authored-By: <your model> <noreply@anthropic.com>"
```

---

### Task 6: Python bindings (`absrbrs.LineList`, `absrbrs.partition`)

**Files:**
- Create: `absrbrs-py/src/lines.rs`, `absrbrs-py/src/partition.rs`, `python/absrbrs/lines.py`, `python/absrbrs/partition.py`
- Modify: `absrbrs-py/src/lib.rs`, `python/absrbrs/__init__.py`, `python/absrbrs/_absrbrs/__init__.pyi`, `tests/test_package.py`, `docs/api.md`
- Test: `tests/test_lines.py`, `tests/test_partition.py`

**Interfaces:**
- Consumes: `absrbrs::lines::{LineList, LineColumns, LineField, molecule_id, molecule_name, isotopologue}`, `absrbrs::io::{read_hitran_par, read_hapi}`, `absrbrs::partition::{tips2021, PartitionFunction}`, and `crate::errors::to_pyerr`.
- Produces in the extension:
  - `class LineList` (`#[pyclass(frozen, module = "absrbrs._absrbrs")]`) with:
    - static methods `from_hitran(path: str)` and `from_hapi(directory: str, table: str)`;
    - `from_columns(molecule: NDArray[uint16], isotopologue: NDArray[uint16], nu, s_ref, e_lower, gamma_air, gamma_self, n_air, delta_air: NDArray[float64], n_self, delta_self, deltap_air: NDArray[float64] | None)`;
    - `__len__`;
    - `column(name: str) -> NDArray[float64] | None`, `molecule() -> NDArray[uint16]` and `isotopologue() -> NDArray[uint16]`, each a fresh copy;
    - `select(molecule: int | None, isotopologues: NDArray[uint16] | None, wavenumber: tuple[float, float] | None) -> LineList`;
    - `with_column(name: str, values: NDArray[float64]) -> LineList`.
  - `tips2021(molecule: int, isotopologue: int, temperature: NDArray[float64]) -> tuple[NDArray[float64], NDArray[float64]]`.
  - `molecule_id(name: str) -> int`, raising InvalidInputError for an unknown name.
  - `isotopologue_info(molecule: int, isotopologue: int) -> dict`, with keys `formula`, `molecule_name`, `abundance`, `mass` and `global_id`; it raises InvalidInputError if the isotopologue is unknown.
- Produces in Python:
  - `absrbrs.LineList`, a wrapper class with numpydoc docstrings. Its methods are:
    - `LineList.from_hitran(path)`, `LineList.from_hapi(directory, table)`;
    - `LineList.from_arrays(*, molecule, isotopologue, nu, s_ref, e_lower, gamma_air, gamma_self, n_air, delta_air, n_self=None, delta_self=None, deltap_air=None)`;
    - `len(lines)`;
    - read-only array properties for every column: `molecule`, `isotopologue`, `nu`, `s_ref`, `e_lower`, `gamma_air`, `gamma_self`, `n_air`, `delta_air`, and `n_self`, `delta_self`, `deltap_air` (None if absent);
    - `select(*, molecule=None, isotopologues=None, wavenumber=None)`, where `molecule` may be an int or a name;
    - `with_column(name, values)`;
    - `__repr__`.
  - `absrbrs.partition.tips2021(molecule, isotopologue, temperature) -> (q, dlnq_dt)`, which broadcasts over `temperature`; `molecule` may be an int or a name.

Binding conventions are as in phase 1a (`absrbrs-py/src/special.rs` is the reference):
- `Bound<PyArray1<T>>` with `try_readonly()`, and BufferError on a failed borrow;
- outputs reserved with `try_reserve_exact`, and MemoryError on failure;
- the GIL released (`py.detach`) for reading files and for evaluating tables;
- Python converts types and shapes, and Rust validates values;
- real-only arguments use dtype kinds `"iuf"`, with TypeError otherwise;
- returned arrays are read-only: in Python, set `array.flags.writeable = False` on each column copy.

- [ ] **Step 1: Write the failing tests**

Create `tests/test_lines.py`:
```python
import numpy as np
import pytest

import absrbrs
from absrbrs import LineList


def make(n=3, **extra):
    kw = dict(
        molecule=[2] * n, isotopologue=list(range(1, n + 1)), nu=[3.0, 1.0, 2.0][:n],
        s_ref=[1e-21, 2e-21, 3e-21][:n], e_lower=[10.0] * n, gamma_air=[0.07] * n,
        gamma_self=[0.09] * n, n_air=[0.75] * n, delta_air=[-0.004] * n,
    )
    kw.update(extra)
    return LineList.from_arrays(**kw)


def test_from_arrays_sorts_and_exposes_read_only_columns():
    lines = make()
    assert len(lines) == 3
    np.testing.assert_array_equal(lines.nu, [1.0, 2.0, 3.0])
    np.testing.assert_array_equal(lines.isotopologue, [2, 3, 1])
    assert lines.isotopologue.dtype == np.uint16
    assert lines.n_self is None
    with pytest.raises(ValueError, match="read-only"):
        lines.nu[0] = 5.0


def test_optional_columns_and_with_column():
    lines = make(n_self=[0.6, np.nan, 0.5])
    np.testing.assert_array_equal(lines.n_self, [np.nan, 0.5, 0.6])  # sorted with the rows
    g = lines.with_column("gamma_air", [0.01, 0.02, 0.03])
    np.testing.assert_array_equal(g.gamma_air, [0.01, 0.02, 0.03])
    np.testing.assert_array_equal(lines.gamma_air, [0.07] * 3)  # the original is unchanged
    with pytest.raises(absrbrs.InvalidInputError):
        lines.with_column("gamma_air", [0.01])
    with pytest.raises(absrbrs.InvalidInputError, match="unknown"):
        lines.with_column("sw", [1, 2, 3])


def test_select_by_molecule_name_isotopologues_and_range():
    lines = make()
    assert len(lines.select(molecule="CO2", isotopologues=[1, 3], wavenumber=(1.5, 3.0))) == 2
    assert len(lines.select(molecule=1)) == 0
    with pytest.raises(absrbrs.InvalidInputError, match="molecule"):
        lines.select(molecule="unobtainium")


def test_invalid_arrays():
    with pytest.raises(absrbrs.InvalidInputError):
        make(gamma_air=[0.07, np.nan, 0.07])
    with pytest.raises(TypeError):
        make(nu=[1 + 1j, 2.0, 3.0])


def fixed(v, width, prec):
    """HITRAN's Fortran fixed format: drop the leading zero if needed (".0700", "-.004000")."""
    s = f"{v:.{prec}f}"
    if len(s) > width:
        s = s.replace("0.", ".", 1)
    return s.rjust(width)


def par_record(mol, iso, nu, s):
    r = (
        f"{mol:2d}{iso}{fixed(nu, 12, 6)}{s:10.3E}{1.0:10.3E}{fixed(0.07, 5, 4)}{fixed(0.09, 5, 3)}"
        f"{fixed(100.0, 10, 4)}{fixed(0.75, 4, 2)}{fixed(-0.004, 8, 6)}"
    )
    assert len(r) == 67
    return r.ljust(160)


def test_from_hitran(tmp_path):
    path = tmp_path / "CO2.par"
    path.write_text(par_record(2, "1", 6300.5, 1e-23) + "\n" + par_record(2, "A", 6200.0, 1e-25) + "\n")
    lines = LineList.from_hitran(path)
    np.testing.assert_array_equal(lines.isotopologue, [11, 1])
    with pytest.raises(FileNotFoundError):
        LineList.from_hitran(tmp_path / "missing.par")
    (tmp_path / "bad.par").write_text("short\n")
    with pytest.raises(absrbrs.ParseError) as info:
        LineList.from_hitran(tmp_path / "bad.par")
    assert info.value.line_number == 1
```

Create `tests/test_partition.py`:
```python
import numpy as np
import pytest

import absrbrs
from absrbrs import partition


def test_tips2021_values_and_derivative():
    q, d = partition.tips2021("H2O", 1, [296.0, 1.0])
    np.testing.assert_allclose(q, [174.58130, 1.0], rtol=1e-7)
    t = np.array([150.25, 296.0, 777.7])
    _, d = partition.tips2021(2, 1, t)
    qp, _ = partition.tips2021(2, 1, t + 1e-4)
    qm, _ = partition.tips2021(2, 1, t - 1e-4)
    np.testing.assert_allclose(d, (np.log(qp) - np.log(qm)) / 2e-4, rtol=1e-6)


def test_shapes():
    q, d = partition.tips2021(2, 1, np.full((2, 3), 296.0))
    assert q.shape == d.shape == (2, 3)
    q, _ = partition.tips2021(2, 1, 296.0)
    assert q.shape == ()


def test_errors():
    with pytest.raises(absrbrs.OutsideTableError) as info:
        partition.tips2021(3, 1, [296.0, 1500.0])
    assert info.value.temperature == 1500.0
    with pytest.raises(absrbrs.MissingPartitionFunctionError):
        partition.tips2021(10, 3, 296.0)
    with pytest.raises(absrbrs.InvalidInputError, match="molecule"):
        partition.tips2021("XYZ", 1, 296.0)
```
In `tests/test_package.py`, extend the expected `_absrbrs.__all__` with `"LineList"`, `"isotopologue_info"`, `"molecule_id"` and `"tips2021"`.

Run: `uv run --locked pytest -q tests/test_lines.py tests/test_partition.py`. Expected: ImportError.

- [ ] **Step 2: Implement the extension and the Python layer**

Implement `absrbrs-py/src/lines.rs` and `absrbrs-py/src/partition.rs` to the interfaces above:
- The `LineList` pyclass wraps `absrbrs::lines::LineList`.
- `column(name)` maps an unknown name to InvalidInputError (`unknown line field "sw"`).
- `tips2021` evaluates every temperature with the GIL released. On the first error it stops and raises the mapped core error.

Register everything in `absrbrs-py/src/lib.rs`.

Then write `python/absrbrs/lines.py`, with the `LineList` wrapper holding `_inner`. Requirements:
- `molecule` is resolved through `_absrbrs.molecule_id` when it is a `str`;
- arrays are converted with `np.asarray(..., dtype=float64)` after a real-only dtype check, and `uint16` for molecule and isotopologue;
- properties return read-only arrays.

Write `python/absrbrs/partition.py`: broadcast `temperature` and flatten it to contiguous float64, resolve the molecule, call the extension, and reshape.

Export `LineList` and `partition` from `absrbrs/__init__.py`, add them to `__all__`, and update the stubs. In `docs/api.md`, add sections for `absrbrs.LineList` (`.. autoclass:: absrbrs.LineList` with `:members:`) and `absrbrs.partition`.

- [ ] **Step 3: Run every Python check and commit**

Run:
```bash
uv sync
uv run --locked pytest -q
uv run --locked python -m mypy.stubtest absrbrs._absrbrs
uv run --locked mypy --strict python/absrbrs/lines.py python/absrbrs/partition.py
uv run --locked ruff check . && uv run --locked ruff format --check .
cargo clippy --workspace --all-targets -- -D warnings && cargo +1.99.0 clippy --workspace --all-targets -- -D warnings
uv run --locked --extra docs sphinx-build -E -W --keep-going -b html docs docs/_build/html
```
Expected: everything passes.

Commit:
```bash
git add absrbrs-py python tests docs/api.md
git commit -m "feat(python): absrbrs.LineList and absrbrs.partition.tips2021

Co-Authored-By: <your model> <noreply@anthropic.com>"
```

---

### Task 7: Cross-check against HAPI, docs and wrap-up

**Files:**
- Modify: `validation/pyproject.toml` (add `hitran-api==1.3.0.0`), `validation/uv.lock`
- Create: `validation/tests/test_tips.py`, `docs/theory/line-state.md`, `docs/theory/partition-functions.md`, `docs/plans/phase-1b-followups.md`
- Modify: `docs/index.md` (Theory toctree), `docs/references.md`, `CHANGELOG.md`, `README.md`, `.github/dependabot.yml` (ignore `hitran-api` in `/validation`), `CONTRIBUTING.md`

**Interfaces:**
- Consumes: `absrbrs.partition.tips2021` and HAPI's `partitionSum(M, I, T, version=2021)`.

- [ ] **Step 1: Compare with HAPI**

Add `"hitran-api==1.3.0.0"` to `validation/pyproject.toml` dependencies; it is pinned exactly, and its dependabot entry ignores it. Then run `uv lock --directory validation`.

Create `validation/tests/test_tips.py`:
```python
"""absrbrs's TIPS-2021 against HAPI's embedded TIPS-2021 (version=2021) at HAPI's own nodes.

HAPI embeds the tables at 10 K steps (1, 10, 20, … K), where its interpolation returns the node
values, so the comparison is of node values: equal to within f32 rounding (absrbrs stores f32)
and the tables' printed digits.
"""

import contextlib
import io

import numpy as np
import pytest

from absrbrs import partition

with contextlib.redirect_stdout(io.StringIO()):
    import hapi


def isotopologues():
    tips = {}
    for m, i in sorted(hapi.ISO):
        try:
            partition.tips2021(m, i, 296.0)
        except Exception:
            continue
        tips[(m, i)] = True
    return sorted(tips)


@pytest.mark.parametrize(("molecule", "isotopologue"), isotopologues())
def test_tips_matches_hapi_at_its_nodes(molecule, isotopologue):
    t = np.array([10.0 * k for k in range(1, 101)])  # 10 … 1000 K, inside every table
    ours, _ = partition.tips2021(molecule, isotopologue, t)
    with contextlib.redirect_stdout(io.StringIO()):
        theirs = np.array([hapi.partitionSum(molecule, isotopologue, float(x), version=2021) for x in t])
    np.testing.assert_allclose(ours, theirs, rtol=2e-7)
```
Run: `uv run --locked --directory validation python -m pytest -q tests/test_tips.py`.

Expected: about 146 tests pass. If a molecule fails at the 2e-7 level, inspect what HAPI's embedded table holds at that node, using `hapi.TIPS_2021_ISOQ_HASH`. Report the case and do not widen the tolerance; for example, HAPI may store fewer digits.

- [ ] **Step 2: Documentation**

Create `docs/theory/line-state.md`. It states the formulas of Task 5 exactly as implemented, and covers:
- the constants (exact SI, with c₂ to the nearest double, and how it compares with sasktran2's 1.4387769);
- the units;
- the optional-field defaults and the tie rule from the Global Constraints;
- the fact that the state and line Jacobians are separate functions, so state-only work never computes line terms.

Create `docs/theory/partition-functions.md`. It covers:
- the TIPS-2021 source and licence;
- the 1 K storage as f32 (≤ 6e-8);
- the cubic Hermite in ln Q, with its consistent derivative;
- each table's range, and `OutsideTable` outside it;
- user tables;
- why HAPI's 10 K tables are not the reference. Phase 1b planning measured them at 3–5e-5 at atmospheric temperatures, and up to 50% below 20 K, against the 1 K source.

Add both pages to the Theory toctree in `docs/index.md`. Add to `docs/references.md`: Kochanov et al. 2016 (HAPI) if it is not already there, and the Zenodo TIPS-2021 record.

- [ ] **Step 3: Follow-ups, changelog, README**

Create `docs/plans/phase-1b-followups.md`. Its sections:
- **Decisions:**
  - TIPS-2021 is stored at 1 K steps as f32. This supersedes spec §5.5's "equal HAPI's partitionSum exactly at nodes": the check is now against the source to f32 rounding, and against HAPI at its nodes to 2e-7.
  - The tie rule for absent optional fields refines spec §5.3.
  - The constants are exact SI.
  - The line-state layer stays crate-private until the phase 1c engine.
- **Measured:** the local reader timings from Task 3, the TIPS package size, the largest error against Arb, and the HAPI comparison.
- **For phase 1c:** the engine resolves masses with `lines::isotopologue` and partition values with `PartitionFunction`. It computes `PartitionRatio` once per (condition, isotopologue).
- **For phase 2:** the AER line lists (local copy in sasktran2's database), and line mixing.

In `CHANGELOG.md`, add entries under `[Unreleased]` → `### Added`. Update the README status to: phase 1b complete; phase 1c (cross-section engine) is next. In `CONTRIBUTING.md`, add the opt-in `ABSRBRS_HITRAN_LINES_DIR` test and the two `tools/` scripts with their regeneration commands.

- [ ] **Step 4: Run every check and commit**

Run the whole "Checks" block of `CONTRIBUTING.md`, plus `cargo +1.99.0 clippy --workspace --all-targets -- -D warnings`, `uv run --locked --directory validation python -m pytest -q` and `uv run --locked pre-commit run --all-files`. Expected: everything passes.

Commit:
```bash
git add validation docs CHANGELOG.md README.md .github/dependabot.yml CONTRIBUTING.md
git commit -m "docs, validation: TIPS-2021 against HAPI; line-state and partition-function theory; phase 1b follow-ups

Co-Authored-By: <your model> <noreply@anthropic.com>"
```

---

## Phase 1b exit checklist (spec §13)

- [ ] `LineList`, the HITRAN and HAPI readers, and the isotopologue table: Tasks 1–3.
- [ ] TIPS-2021 partition functions with consistent d ln Q/dT, plus user tables: Task 4.
- [ ] Line-state layer with closed-form ∂q/∂θ, checked by dual numbers and against Arb: Task 5.
- [ ] Python: `absrbrs.LineList`, `absrbrs.partition`: Task 6.
- [ ] Licence conditions: the TIPS-2021 and HAPI notices are in `license-files`; no HITRAN data is committed. Tasks 1, 4.
- [ ] Docs, cross-check against HAPI, and follow-ups: Task 7.
