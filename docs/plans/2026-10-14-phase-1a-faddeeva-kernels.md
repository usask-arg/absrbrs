# absrbrs Phase 1a (Faddeeva Kernels) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use usask-foundation:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement the Faddeeva function w(z) and its derivative in Rust at both accuracy tiers, together with the Voigt kernel triplet (K, K_x, K_y), and expose them to Python. Validate them against Arb in CI and benchmark them against sasktran2's J-pole approximation. This is phase 1a of `docs/specs/2026-10-14-absrbrs-design.md` (§6, §13).

**Architecture:** A new `absrbrs::faddeeva` module provides a small public API: `Tier`, `w`, `Kernel`, `kernel` and `kernel_into`. Three private submodules sit behind it:

- `accurate`: a port of Steven G. Johnson's MIT-licensed Faddeeva Package. It provides `w(z)`, `erfcx(x)` and `w_im(x)` (Dawson's integral times 2/√π). The derivative uses the three-way rule measured in phase 1a planning.
- `fast`: the composite chosen in phase 0 (`docs/plans/phase-0-followups.md`). It uses continued-fraction bands, exp(−z²) plus a continued fraction near the real axis, a small-y Taylor series and a Weideman-24 core.
- `cf`: the Laplace continued fraction with its derivative, shared by both tiers.

Both tiers reflect x < 0 onto x > 0, so the symmetries hold exactly. Python gets `absrbrs.special`. A new `validation/` uv project computes dense Arb maps when CI runs and checks both tiers against them through the Python bindings.

**Tech Stack:**
- Rust (stable, edition 2024, MSRV 1.85), with num-complex 0.4. Dev-dependencies: proptest 1, criterion 0.7.
- pyo3 and numpy 0.29.
- python-flint 0.9.0 (Arb), uv, Sphinx.

## Global Constraints

Each item is quoted from the spec and the phase 0 follow-ups, and every task implicitly includes them.

- **Core crate API rules.** No pyo3 or ndarray in the core public API. `#![forbid(unsafe_code)]` and `#![warn(missing_docs)]`. No panics on finite input. One `thiserror` error enum, marked `#[non_exhaustive]`. A new variant also needs a mapping in `absrbrs-py/src/errors.rs`, `python/absrbrs/errors.py` and the `_testing` hook.
- **Dependencies.** Every new dependency needs the usask-foundation:dependency-audit skill. Phase 1a adds:
  - `num-complex` 0.4 (core);
  - `proptest` 1 and `criterion` 0.7 (dev);
  - `numpy` 0.29 (bindings).

  scatrs at commit 4793160 already uses these exact versions; cite its use in each audit.
- **Validated domain (spec §6.2).**
  - y ∈ [0, 10⁶], with y = 0 included, and |x| ≤ 10⁸.
  - Value error is relative to K = Re w. Points where the true K is below `f64::MIN_POSITIVE` are exempt, unless the result is non-finite.
  - Derivative error is componentwise, max(|ΔRe w′|, |ΔIm w′|), relative to |w′|.
- **Tiers (spec §6.3):**

  | Tier | Value error | Derivative error |
  |---|---|---|
  | Fast (the default) | ≤ 1 × 10⁻⁶ | ≤ 1 × 10⁻⁵ |
  | Accurate | ≤ 2 × 10⁻¹³ | ≤ 2 × 10⁻¹² |

- **Kernel triplet.** K = Re w, K_x = Re w′, K_y = −Im w′ (spec §6.1). In the Fast tier's Taylor region, w′ stands for K_x − iK_y of the formula itself.
- **Faddeeva Package licence conditions** (`THIRD_PARTY.md`; row confirmed 2026-10-14):
  - Port from `Faddeeva.cc` only, revision of 12 May 2015, with SHA-256 `46dc0d4e7d81ddcce5a95953b97d31696ffb424228e751aac10dfd9aeb3430a0`, downloaded from `http://ab-initio.mit.edu/Faddeeva.cc` (the site refuses HTTPS).
  - Never port from the ACM code of TOMS Algorithms 916 or 680.
  - Every ported file carries the "Copyright (c) 2012 Massachusetts Institute of Technology" notice and the MIT permission notice in its header.
  - In the same commit as the port, add `absrbrs/LICENSES/Faddeeva-MIT.txt` and list it in `pyproject.toml` `license-files`.
- **Weideman (1994):** derive the coefficients from the paper's equations; do not copy its code listing.
- **Copying from scatrs:** only with `git -C ~/dev/scat-rs show 4793160:<path>`. Never read from or write to the scat-rs working tree.
- **Formatting before checks.** After writing Rust, run `cargo fmt --all`; after writing Python, run `uv run --locked ruff format .`. The plan's code blocks are not always formatted, and the checks run `--check`.
- **Commits** end with a `Co-Authored-By:` trailer naming the model that wrote the commit.
- **Remote:** do not push or add a remote.

## File Map

```
Cargo.toml                                   workspace deps: num-complex            Task 1
absrbrs/Cargo.toml                           deps, dev-deps, bench, exclude         Tasks 1, 3, 4
absrbrs/LICENSES/Faddeeva-MIT.txt            notice for the port                    Task 1
absrbrs/src/lib.rs                           pub mod faddeeva; mod test_data        Task 1
absrbrs/src/test_data.rs                     reference files for unit tests         Task 1
absrbrs/src/faddeeva/mod.rs                  public API                             Tasks 1, 2, 3
absrbrs/src/faddeeva/accurate/{mod,erfcx,w_im}.rs   port: erfcx, w_im                Task 1
absrbrs/src/faddeeva/accurate/w.rs           port: w(z)                             Task 2
absrbrs/src/faddeeva/cf.rs                   continued fraction (shared)            Task 2
absrbrs/src/faddeeva/fast/{mod,weideman,taylor}.rs   Fast tier                       Task 3
absrbrs/tests/faddeeva_accuracy.rs           both tiers vs committed Arb grid       Tasks 2, 3
absrbrs/tests/faddeeva_fast_vs_accurate.rs   proptest + boundary sweep              Task 3
absrbrs/benches/faddeeva.rs                  criterion                              Task 4
reference/absrbrs_reference/generators/weideman.py, reference/tests/test_weideman.py, data   Task 3
absrbrs-py/{Cargo.toml, src/lib.rs, src/special.rs}                                 Task 5
python/absrbrs/{special.py, __init__.py, _absrbrs/__init__.pyi}, tests/test_special.py      Task 5
validation/ (uv project), .github/workflows/ci.yml, .github/dependabot.yml          Task 6
docs/{theory/faddeeva.md, validation/faddeeva-accuracy.md, index.md, api.md}        Tasks 5, 6
docs/plans/phase-1a-followups.md, CHANGELOG.md, pyproject.toml                      Tasks 1, 4, 6
.github/workflows/{ci,scheduled}.yml (ABSRBRS_REQUIRE_REFERENCE_DATA)               Task 1
```

Every command runs from the repository root unless a step says otherwise. Work on branch `phase-1a-faddeeva`; it already exists and is checked out.

---

### Task 1: Accurate tier, real-axis functions (`erfcx`, `w_im`)

**Files:**
- Modify: `Cargo.toml`, adding `num-complex = "0.4.6"` to `[workspace.dependencies]`
- Modify: `absrbrs/Cargo.toml`, adding `num-complex.workspace = true` to `[dependencies]`
- Create: `absrbrs/LICENSES/Faddeeva-MIT.txt`
- Modify: `pyproject.toml:12`, `license-files`
- Create: `absrbrs/src/test_data.rs` (copied from scatrs)
- Modify: `absrbrs/src/lib.rs`
- Create: `absrbrs/src/faddeeva/mod.rs`, `absrbrs/src/faddeeva/accurate/mod.rs`, `absrbrs/src/faddeeva/accurate/erfcx.rs`, `absrbrs/src/faddeeva/accurate/w_im.rs`
- Modify: `.github/workflows/ci.yml` and `.github/workflows/scheduled.yml` (reference-data environment variable)

**Interfaces:**
- Produces:
  - `pub(crate) fn crate::faddeeva::accurate::erfcx(x: f64) -> f64`, which is exp(x²) erfc(x) for all real x;
  - `pub(crate) fn crate::faddeeva::accurate::w_im(x: f64) -> f64`, which is Im w(x) = (2/√π) D(x) for real x, D being Dawson's integral;
  - `pub(crate) fn crate::test_data::reference_file(name: &str) -> Option<String>`;
  - `pub(crate) const FRAC_1_SQRT_PI: f64` in `crate::faddeeva`.

- [ ] **Step 1: Audit the dependency and get the source**

Run the usask-foundation:dependency-audit skill for `num-complex` 0.4. It is MIT OR Apache-2.0, part of rust-num, and scatrs already depends on it with no unsafe-code concerns for our use. Record the outcome in the commit message body.

Download the source and check its hash:
```bash
mkdir -p /tmp/absrbrs-port && curl -s -A "Mozilla/5.0" -o /tmp/absrbrs-port/Faddeeva.cc http://ab-initio.mit.edu/Faddeeva.cc
shasum -a 256 /tmp/absrbrs-port/Faddeeva.cc
```
Expected: `46dc0d4e7d81ddcce5a95953b97d31696ffb424228e751aac10dfd9aeb3430a0`. On any other hash, stop and report.

- [ ] **Step 2: Write the licence notice and wire it in**

Create `absrbrs/LICENSES/Faddeeva-MIT.txt` with these parts, in order:
1. The heading `Faddeeva Package (absrbrs/src/faddeeva/accurate/ is a Rust port of Faddeeva.cc)`.
2. A blank line.
3. The copyright and permission notice exactly as it appears in the header of `Faddeeva.cc` (lines 1–21 of the file: the "Copyright (c) 2012 Massachusetts Institute of Technology" line and the MIT permission text). Copy the text without the `/*`, `*/` and leading ` *` comment markers.

Set `pyproject.toml` line 12 to:
```toml
license-files = ["LICENSE", "absrbrs/LICENSES/Faddeeva-MIT.txt"]
```

- [ ] **Step 3: Copy the reference-file helper and restore the CI variable**

```bash
git -C ~/dev/scat-rs show 4793160:scatrs/src/test_data.rs | sed 's/SCATRS_REQUIRE_REFERENCE_DATA/ABSRBRS_REQUIRE_REFERENCE_DATA/g; s/scatrs/absrbrs/g' > absrbrs/src/test_data.rs
```
`reference_file(name)` returns the text of `reference/data/<name>`, or `None` when the file is missing, as it is in the published crate. With `ABSRBRS_REQUIRE_REFERENCE_DATA=1` a missing file panics instead.

Add the variable back to CI:
- In `.github/workflows/ci.yml` and `.github/workflows/scheduled.yml`, insert after the top-level `permissions:` block:
  ```yaml
  env:
    # Unit tests that read reference/data skip without it (as in the published crate);
    # here a missing file is an error.
    ABSRBRS_REQUIRE_REFERENCE_DATA: "1"
  ```
- In `ci.yml`'s step `Test the packaged crate on its own`, add, before `run: |`:
  ```yaml
          env:
            ABSRBRS_REQUIRE_REFERENCE_DATA: "0"
  ```

- [ ] **Step 4: Write the module skeleton and the failing tests**

Replace `absrbrs/src/lib.rs` with:
```rust
//! Molecular absorption cross sections with first-class linearization.
//!
//! The design is described in `docs/specs/2026-10-14-absrbrs-design.md` in the repository.
#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod faddeeva;
#[cfg(test)]
mod test_data;

pub use error::{Error, Result};
```

Create `absrbrs/src/faddeeva/mod.rs`. The `allow`s are temporary: Task 2 replaces this file with
the public API, which uses everything.
```rust
//! The Faddeeva function w(z) = exp(−z²) erfc(−iz) and the Voigt kernel (spec §6).

#[allow(dead_code)] // used by the public API from Task 2 on
mod accurate;

/// 1/√π.
#[allow(dead_code)] // used from Task 2 on
pub(crate) const FRAC_1_SQRT_PI: f64 = 0.564_189_583_547_756_3;
```

Create `absrbrs/src/faddeeva/accurate/mod.rs`:
```rust
//! The Accurate tier: a Rust port of Steven G. Johnson's Faddeeva Package (`Faddeeva.cc`,
//! revision of 12 May 2015), MIT licence; see `absrbrs/LICENSES/Faddeeva-MIT.txt`.

mod erfcx;
mod w_im;

pub(crate) use erfcx::erfcx;
pub(crate) use w_im::w_im;
```

Create `absrbrs/src/faddeeva/accurate/erfcx.rs` containing only the header and tests for now. The header is the comment below, followed by the full MIT notice from Step 2 as `//` comment lines:
```rust
//! The scaled complementary error function erfcx(x) = exp(x²) erfc(x), for real x.
//!
//! Ported from `Faddeeva.cc` (revision of 12 May 2015): `erfcx_y100` and `FADDEEVA_RE(erfcx)`.
// Copyright (c) 2012 Massachusetts Institute of Technology
// <the MIT permission notice from Step 2, as // comment lines>
#![allow(clippy::excessive_precision)] // the coefficients are copied verbatim

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_data::reference_file;

    #[test]
    fn known_values() {
        assert_eq!(erfcx(0.0), 1.0);
        assert!((erfcx(1.0) - 0.427_583_576_155_807).abs() <= 1e-16);
    }

    /// The x = 0 column of the Faddeeva grid: w(iy) = erfcx(y), for y from 0 to 1e6.
    #[test]
    fn matches_arb_on_the_imaginary_axis() {
        let Some(text) = reference_file("faddeeva.jsonl") else { return };
        let mut worst = 0.0f64;
        for line in text.lines().filter(|l| !l.trim().is_empty()) {
            let r: serde_json::Value = serde_json::from_str(line).unwrap();
            let y = r["inputs"]["y"].as_f64().unwrap();
            let expected = r["outputs"]["w_re"]["hi"][0].as_f64().unwrap(); // x[0] = 0
            let err = ((erfcx(y) - expected) / expected).abs();
            worst = if err.is_nan() { f64::INFINITY } else { worst.max(err) };
        }
        eprintln!("erfcx: worst relative error {worst:e}");
        assert!(worst <= 1e-14, "{worst:e}");
    }

    /// erfcx(−x) = 2 exp(x²) − erfcx(x).
    #[test]
    fn negative_arguments() {
        for x in [0.25, 0.5, 1.0, 3.0, 10.0, 26.0] {
            let expected = 2.0 * (x * x).exp() - erfcx(x);
            assert!(((erfcx(-x) - expected) / expected).abs() <= 1e-14, "x = {x}");
        }
    }
}
```

Create `absrbrs/src/faddeeva/accurate/w_im.rs` the same way: the header (ported functions `w_im_y100` and `FADDEEVA(w_im)`, the same notice and the same `allow`), followed by:
```rust
#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_data::reference_file;

    #[test]
    fn known_values() {
        assert_eq!(w_im(0.0), 0.0);
        assert!((w_im(1.0) - 0.607_157_705_841_393_7).abs() <= 1e-16);
    }

    /// The y = 0 record of the Faddeeva grid: Im w(x) = w_im(x), for x from 0 to 1e8.
    #[test]
    fn matches_arb_on_the_real_axis() {
        let Some(text) = reference_file("faddeeva.jsonl") else { return };
        let first = text.lines().next().unwrap(); // the record for y = 0
        let r: serde_json::Value = serde_json::from_str(first).unwrap();
        assert_eq!(r["inputs"]["y"].as_f64(), Some(0.0));
        let xs = r["inputs"]["x"].as_array().unwrap();
        let mut worst = 0.0f64;
        for (i, x) in xs.iter().enumerate() {
            let x = x.as_f64().unwrap();
            let expected = r["outputs"]["w_im"]["hi"][i].as_f64().unwrap();
            let got = w_im(x);
            let err = if expected == 0.0 { got.abs() } else { ((got - expected) / expected).abs() };
            worst = if err.is_nan() { f64::INFINITY } else { worst.max(err) };
        }
        eprintln!("w_im: worst relative error {worst:e}");
        assert!(worst <= 1e-14, "{worst:e}");
    }

    #[test]
    fn is_odd() {
        for x in [1e-8, 0.3, 1.0, 5.5, 45.0, 1e5] {
            assert_eq!(w_im(-x), -w_im(x));
        }
    }
}
```

Run: `cargo test -p absrbrs faddeeva`

Expected: compile errors ``cannot find function `erfcx` `` and ``cannot find function `w_im` ``.

- [ ] **Step 5: Port the functions**

Translate the following from `/tmp/absrbrs-port/Faddeeva.cc` into the two files, above the test modules:
- **Into `erfcx.rs`:**
  - `static double erfcx_y100(double y100)` (lines 1013–1421) as `fn erfcx_y100(y100: f64) -> f64`;
  - `double FADDEEVA_RE(erfcx)(double x)` (lines 1422–1457) as `pub(crate) fn erfcx(x: f64) -> f64`.
- **Into `w_im.rs`:**
  - `static double w_im_y100(double y100, double x)` (lines 1458–1863) as `fn w_im_y100(y100: f64, x: f64) -> f64`;
  - `double FADDEEVA(w_im)(double x)` (lines 1864–1900) as `pub(crate) fn w_im(x: f64) -> f64`.

Port rules, which make the translation line-by-line reviewable:
- Keep the C `switch` statements as `match` on the integer the C code switches on, in the same order and with the same arms.
- Copy every numeric literal verbatim, including digits beyond f64 precision.
- Keep the comments that explain the algorithm.
- Map the library calls and constants as follows:

  | C | Rust |
  |---|---|
  | `exp` | `f64::exp` |
  | `floor` | `f64::floor` |
  | `copysign` | `f64::copysign` |
  | `isnan` | `is_nan` |
  | `Inf` | `f64::INFINITY` |
  | `NaN` | `f64::NAN` |

  `ispi` becomes a local `const` with Johnson's literal `0.56418958354775628694807945156`.
- Integer conversions `(int) y100` become `y100 as i32`. y100 is bounded by the surrounding range checks, as in the C code.
- Do not "improve" anything. Any deviation from the C code needs a comment saying why.

- [ ] **Step 6: Run the tests and lint checks**

Run:
```bash
cargo test -p absrbrs faddeeva -- --nocapture
cargo clippy --workspace --all-targets -- -D warnings
cargo fmt --all -- --check
```
Expected:
- `6 passed` for the filter `faddeeva`. The printed worst errors should be around 1e-16 to 1e-15.
- clippy and fmt clean.

If an error exceeds 1e-14, compare the failing range's code with the C source line by line. Do not loosen the tolerance.


- [ ] **Step 7: Commit**

```bash
git add Cargo.toml Cargo.lock absrbrs/Cargo.toml absrbrs/LICENSES absrbrs/src pyproject.toml .github/workflows/ci.yml .github/workflows/scheduled.yml
git commit -m "feat(faddeeva): port erfcx and w_im from the Faddeeva Package (Accurate tier)

Dependency audit: num-complex 0.4 (MIT OR Apache-2.0) ...

Co-Authored-By: <your model> <noreply@anthropic.com>"
```

---

### Task 2: Accurate tier, w(z) with its derivative, and the public API

**Files:**
- Create: `absrbrs/src/faddeeva/accurate/w.rs`, `absrbrs/src/faddeeva/cf.rs`
- Modify: `absrbrs/src/faddeeva/accurate/mod.rs`, `absrbrs/src/faddeeva/mod.rs`
- Test: `absrbrs/tests/faddeeva_accuracy.rs`

**Interfaces:**
- Consumes: `accurate::erfcx`, `accurate::w_im` and `FRAC_1_SQRT_PI` from Task 1.
- Produces (public, in `absrbrs::faddeeva`):
  - `pub enum Tier { Accurate }`, with `#[non_exhaustive]` and derives `Clone, Copy, Debug, PartialEq, Eq, Hash`. Task 3 adds `Fast` and makes it the default.
  - `pub fn w(z: Complex64, tier: Tier) -> (Complex64, Complex64)`, returning w(z) and w′(z). Both are NaN when `!(z.im >= 0.0)`.
  - `pub struct Kernel { pub k: f64, pub k_x: f64, pub k_y: f64 }`, deriving `Clone, Copy, Debug, Default, PartialEq`.
  - `pub fn kernel(x: f64, y: f64, tier: Tier) -> Kernel`.
  - `pub fn kernel_into(x: &[f64], y: f64, tier: Tier, k: &mut [f64], k_x: &mut [f64], k_y: &mut [f64]) -> crate::Result<()>`.
  - The re-export `pub use num_complex::Complex64`.
- Produces (crate-private):
  - `cf::laplace_cf(z: Complex64, depth: u32) -> (Complex64, Complex64)`;
  - `cf::gauss_cf(z: Complex64, depth: u32) -> (Complex64, Complex64)`;
  - `accurate::w_and_derivative(z: Complex64) -> (Complex64, Complex64)`, for Re z ≥ 0 and Im z ≥ 0.

The Accurate derivative rule was measured on 564 696 Arb points before this plan was written: the committed grid, the prototype's dense map and its holdout. Johnson's w meets the value target there (9.7e-14). The derivative uses three regions:

| Region | Derivative from | Measured error |
|---|---|---|
| \|z\| < 6 | the identity w′ = −2zw + 2i/√π | 9.3e-13 |
| x ≥ 6 and y ≤ 0.1 | derivative of exp(−z²) + continued fraction, depth 16 | 2.8e-13 |
| elsewhere (\|z\| ≥ 6) | derivative of the continued fraction, depth 20 | 1.7e-13 |

The identity alone fails for |z| ≥ 6, reaching 7.4e-12 near the axis and far more in the wings.

- [ ] **Step 1: Write the failing tests**

Create `absrbrs/tests/faddeeva_accuracy.rs`:
```rust
//! The Faddeeva tiers against the committed Arb grid (reference/data/faddeeva.jsonl), with the
//! error metrics of spec §6.2, and the API's symmetry and domain rules.

mod common;

use absrbrs::faddeeva::{self, Complex64, Kernel, Tier};
use common::load;

/// The largest value and derivative errors of `tier` on the grid, with where they occur.
fn worst(tier: Tier) -> ((f64, f64, f64), (f64, f64, f64)) {
    let (mut value, mut derivative) = ((0.0, 0.0, 0.0), (0.0, 0.0, 0.0));
    for r in load("faddeeva.jsonl") {
        let y = r.input_f64("y");
        let (w_re, dw_re, dw_im) = (r.output("w_re"), r.output("dw_re"), r.output("dw_im"));
        for (i, x) in r.input_f64s("x").into_iter().enumerate() {
            let (w, dw) = faddeeva::w(Complex64::new(x, y), tier);
            let k = w_re.hi[i];
            let ve = if k < f64::MIN_POSITIVE {
                if w.re.is_finite() { 0.0 } else { f64::INFINITY }
            } else {
                (w.re - k).abs() / k
            };
            let de = (dw.re - dw_re.hi[i]).abs().max((dw.im - dw_im.hi[i]).abs())
                / dw_re.hi[i].hypot(dw_im.hi[i]);
            let (ve, de) = (nan_is_infinite(ve), nan_is_infinite(de));
            if ve > value.0 {
                value = (ve, x, y);
            }
            if de > derivative.0 {
                derivative = (de, x, y);
            }
        }
    }
    (value, derivative)
}

fn nan_is_infinite(e: f64) -> f64 {
    if e.is_nan() { f64::INFINITY } else { e }
}

#[test]
fn accurate_tier_meets_its_targets() {
    let (v, d) = worst(Tier::Accurate);
    eprintln!("Accurate: value {:e} at ({}, {}), derivative {:e} at ({}, {})", v.0, v.1, v.2, d.0, d.1, d.2);
    assert!(v.0 <= 2e-13, "value error {:e} at x = {}, y = {}", v.0, v.1, v.2);
    assert!(d.0 <= 2e-12, "derivative error {:e} at x = {}, y = {}", d.0, d.1, d.2);
}

#[test]
fn negative_x_follows_the_symmetries_exactly() {
    for tier in [Tier::Accurate] {
        for (x, y) in [(0.3, 0.0), (2.0, 0.5), (7.0, 1e-3), (40.0, 2.0), (1e5, 1e-20)] {
            let (w, dw) = faddeeva::w(Complex64::new(x, y), tier);
            let (wn, dwn) = faddeeva::w(Complex64::new(-x, y), tier);
            assert_eq!(wn, w.conj(), "w at ({x}, {y})");
            assert_eq!(dwn, -dw.conj(), "w′ at ({x}, {y})");
        }
    }
}

#[test]
fn negative_or_nan_y_gives_nan() {
    for z in [Complex64::new(1.0, -1e-300), Complex64::new(1.0, f64::NAN)] {
        let (w, dw) = faddeeva::w(z, Tier::Accurate);
        assert!(w.re.is_nan() && w.im.is_nan() && dw.re.is_nan() && dw.im.is_nan());
    }
}

#[test]
fn kernel_is_the_triplet_of_w() {
    for (x, y) in [(0.0, 0.0), (1.5, 0.25), (-9.0, 1e-3), (300.0, 4.0)] {
        let (w, dw) = faddeeva::w(Complex64::new(x, y), Tier::Accurate);
        assert_eq!(faddeeva::kernel(x, y, Tier::Accurate), Kernel { k: w.re, k_x: dw.re, k_y: -dw.im });
    }
}

#[test]
fn kernel_into_matches_kernel_bitwise() {
    let x: Vec<f64> = (-400..=400).map(|i| f64::from(i) * 0.137).collect();
    for y in [0.0, 1e-3, 0.5, 30.0] {
        let n = x.len();
        let (mut k, mut kx, mut ky) = (vec![0.0; n], vec![0.0; n], vec![0.0; n]);
        faddeeva::kernel_into(&x, y, Tier::Accurate, &mut k, &mut kx, &mut ky).unwrap();
        for i in 0..n {
            let one = faddeeva::kernel(x[i], y, Tier::Accurate);
            assert_eq!((k[i], kx[i], ky[i]), (one.k, one.k_x, one.k_y), "x = {}, y = {y}", x[i]);
        }
    }
}

#[test]
fn kernel_into_rejects_mismatched_lengths() {
    let (mut k, mut kx, mut ky) = (vec![0.0; 3], vec![0.0; 3], vec![0.0; 2]);
    let err = faddeeva::kernel_into(&[0.0, 1.0, 2.0], 1.0, Tier::Accurate, &mut k, &mut kx, &mut ky);
    assert!(matches!(err, Err(absrbrs::Error::InvalidInput(_))));
}
```

Run: `cargo test -p absrbrs --test faddeeva_accuracy`

Expected: a compile error, ``cannot find type `Tier` in module `faddeeva` ``, along with the other missing items.

- [ ] **Step 2: Write the continued fraction**

Create `absrbrs/src/faddeeva/cf.rs`. It is shared by both tiers and has no third-party origin; it follows `laplace_cf` and `gauss_cf` in `prototypes/faddeeva/faddeeva_prototype/candidates.py`.
```rust
//! The Laplace continued fraction for w(z), and exp(−z²) plus it, each with the derivative of
//! its own formula (forward differentiation of the recurrence).

use num_complex::Complex64;

use super::FRAC_1_SQRT_PI;

const I: Complex64 = Complex64::new(0.0, 1.0);

/// w(z) ≈ (i/√π) / (z − ½/(z − 1/(z − (3/2)/(z − …)))) truncated after `depth` levels, and the
/// derivative of that truncation. Depth 0 is the Lorentzian i/(√π z); depth 1 is Humlíček's
/// region-I formula. An asymptotic expansion: for large |z| only (its truncations have poles on
/// the real axis within |x| < 6 for depth ≤ 20).
pub(crate) fn laplace_cf(z: Complex64, depth: u32) -> (Complex64, Complex64) {
    let mut f = z;
    let mut fp = Complex64::new(1.0, 0.0);
    for k in (1..=depth).rev() {
        let a = 0.5 * f64::from(k);
        let inv = f.inv();
        (f, fp) = (z - a * inv, 1.0 + a * fp * inv * inv);
    }
    let inv = f.inv();
    (I * FRAC_1_SQRT_PI * inv, -I * FRAC_1_SQRT_PI * fp * inv * inv)
}

/// exp(−z²) plus the continued fraction, and its derivative: keeps the Gaussian core that the
/// asymptotic expansion lacks. exp(−z²) overflows for y above about 26.6: small y only.
pub(crate) fn gauss_cf(z: Complex64, depth: u32) -> (Complex64, Complex64) {
    let g = (-(z * z)).exp();
    let (w, dw) = laplace_cf(z, depth);
    (w + g, dw - 2.0 * z * g)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shallow_fractions_are_the_classical_formulas() {
        let z = Complex64::new(3.0, 2.0);
        let (w0, dw0) = laplace_cf(z, 0);
        assert!((w0 - I * FRAC_1_SQRT_PI / z).norm() <= 1e-16 * w0.norm());
        assert!((dw0 + I * FRAC_1_SQRT_PI / (z * z)).norm() <= 1e-16 * dw0.norm());
        let (w1, _) = laplace_cf(z, 1);
        let expected = I * FRAC_1_SQRT_PI * z / (z * z - 0.5);
        assert!((w1 - expected).norm() <= 1e-15 * w1.norm());
    }

    /// The derivative is that of the formula: central differences agree.
    #[test]
    fn derivative_is_that_of_the_formula() {
        for (depth, z) in [(20, Complex64::new(9.0, 3.0)), (16, Complex64::new(7.0, 0.01))] {
            let h = 1e-5 * z.norm();
            let fd = (laplace_cf(z + h, depth).0 - laplace_cf(z - h, depth).0) / (2.0 * h);
            let dw = laplace_cf(z, depth).1;
            assert!((fd - dw).norm() <= 1e-8 * dw.norm());
            let fd = (gauss_cf(z + h, depth).0 - gauss_cf(z - h, depth).0) / (2.0 * h);
            let dw = gauss_cf(z, depth).1;
            assert!((fd - dw).norm() <= 1e-8 * dw.norm());
        }
    }
}
```

- [ ] **Step 3: Port w(z) and add the derivative rule**

Create `absrbrs/src/faddeeva/accurate/w.rs` with the same header as Task 1's ported files (function `FADDEEVA(w)` and its helpers). Port from `Faddeeva.cc`:
- the helpers `sinc` (line 609), `sinh_taylor` (614) and `sqr` (619);
- the table `expa2n2` (lines 623–679), as `const EXPA2N2: [f64; 52]`;
- `cmplx FADDEEVA(w)(cmplx z, double relerr)` (lines 680–1012), as `pub(crate) fn w(z: Complex64) -> Complex64`.

Rules on top of Task 1's port rules:
- **Precision path.** Port only the machine-precision path. `relerr` is fixed at `DBL_EPSILON`: keep the branch `if (relerr <= DBL_EPSILON)` with its consequences (`a = 0.518321480430085929872`, `c = 0.329973702884629072537`, `a2 = 0.268657157075235951582`) and the precomputed-table loops. Drop the parameter, the `else` branch that computes `a, c, a2` from `relerr`, and the loops that compute `exp(-a2*(n*n))` on the fly.
- **Continued fraction.** Port only the `#if USE_CONTINUED_FRACTION` branch (the macro is 1). Drop the `#else` branch.
- **Calls.** `FADDEEVA_RE(erfcx)` becomes `super::erfcx`, and `FADDEEVA(w_im)` becomes `super::w_im`.
- **Complex operations.** `cexp(C(a, b))` becomes `Complex64::new(a, b).exp()`, `C(a, b)` becomes `Complex64::new(a, b)`, and `creal`/`cimag` become `.re`/`.im`.
- **The y < 0 branches.** Port them faithfully, so the function can be reviewed against the source. Callers in this crate only pass Re z ≥ 0 and Im z ≥ 0.

Then replace `absrbrs/src/faddeeva/accurate/mod.rs` with:
```rust
//! The Accurate tier: a Rust port of Steven G. Johnson's Faddeeva Package (`Faddeeva.cc`,
//! revision of 12 May 2015), MIT licence; see `absrbrs/LICENSES/Faddeeva-MIT.txt`. The
//! derivative is absrbrs's own (`w_and_derivative`).

mod erfcx;
mod w;
mod w_im;

use num_complex::Complex64;

pub(crate) use erfcx::erfcx;
pub(crate) use w_im::w_im;

use super::FRAC_1_SQRT_PI;
use super::cf::{gauss_cf, laplace_cf};

/// w(z) and w′(z) for Re z ≥ 0 and Im z ≥ 0.
///
/// The value is Johnson's. The derivative uses the identity w′ = −2zw + 2i/√π for |z| < 6,
/// where its cancellation (about |z|² in relative terms) stays within the tier's target. Beyond
/// that it uses the derivative of an asymptotic formula: exp(−z²) plus a depth-16 continued
/// fraction near the real axis (x ≥ 6, y ≤ 0.1), and a depth-20 continued fraction elsewhere.
/// Measured against Arb before implementation: 9.3e-13, 2.8e-13 and 1.7e-13 relative to |w′|
/// in the three regions.
pub(crate) fn w_and_derivative(z: Complex64) -> (Complex64, Complex64) {
    let w = w::w(z);
    let dw = if z.re >= 6.0 && z.im <= 0.1 {
        gauss_cf(z, 16).1
    } else if z.norm() >= 6.0 {
        laplace_cf(z, 20).1
    } else {
        -2.0 * z * w + Complex64::new(0.0, 2.0 * FRAC_1_SQRT_PI)
    };
    (w, dw)
}
```

- [ ] **Step 4: Write the public API**

Replace `absrbrs/src/faddeeva/mod.rs` with:
```rust
//! The Faddeeva function w(z) = exp(−z²) erfc(−iz) and the Voigt kernel (spec §6).
//!
//! Every function takes Im z ≥ 0 (y ≥ 0): the Voigt profile is K(x, y) = Re w(x + iy) with
//! y = γ_L/σ_D ≥ 0. Negative x uses the exact symmetries w(−x + iy) = conj(w(x + iy)) and
//! w′(−x + iy) = −conj(w′(x + iy)). For Im z < 0 or NaN the results are NaN. Accuracy is
//! validated for y ≤ 10⁶ and |x| ≤ 10⁸ (spec §6.2).

mod accurate;
mod cf;

pub use num_complex::Complex64;

use crate::{Error, Result};

/// 1/√π.
pub(crate) const FRAC_1_SQRT_PI: f64 = 0.564_189_583_547_756_3;

/// The accuracy tier of the Faddeeva function (spec §6.3).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum Tier {
    /// Value error ≤ 2e-13 relative to Re w; derivative error ≤ 2e-12 relative to |w′|. A port of
    /// the Faddeeva Package.
    Accurate,
}

/// The Voigt kernel triplet at x + iy: K = Re w, K_x = ∂K/∂x and K_y = ∂K/∂y.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Kernel {
    /// K = Re w(x + iy).
    pub k: f64,
    /// ∂K/∂x = Re w′.
    pub k_x: f64,
    /// ∂K/∂y = −Im w′.
    pub k_y: f64,
}

/// w(z) and its derivative w′(z), for Im z ≥ 0.
pub fn w(z: Complex64, tier: Tier) -> (Complex64, Complex64) {
    if !(z.im >= 0.0) {
        let nan = Complex64::new(f64::NAN, f64::NAN);
        return (nan, nan);
    }
    let reflect = z.re < 0.0;
    let z = if reflect { Complex64::new(-z.re, z.im) } else { z };
    let (w, dw) = match tier {
        Tier::Accurate => accurate::w_and_derivative(z),
    };
    if reflect { (w.conj(), -dw.conj()) } else { (w, dw) }
}

/// The kernel triplet at x + iy, for y ≥ 0.
pub fn kernel(x: f64, y: f64, tier: Tier) -> Kernel {
    let (w, dw) = w(Complex64::new(x, y), tier);
    Kernel { k: w.re, k_x: dw.re, k_y: -dw.im }
}

/// The kernel triplet at every x of `x` for one y, into `k`, `k_x` and `k_y` (each as long as
/// `x`). Results equal those of [`kernel`] bitwise. Returns [`Error::InvalidInput`] if the
/// lengths differ.
pub fn kernel_into(
    x: &[f64],
    y: f64,
    tier: Tier,
    k: &mut [f64],
    k_x: &mut [f64],
    k_y: &mut [f64],
) -> Result<()> {
    let n = x.len();
    if k.len() != n || k_x.len() != n || k_y.len() != n {
        return Err(Error::InvalidInput(format!(
            "kernel_into: x has {n} values but k, k_x and k_y have {}, {} and {}",
            k.len(),
            k_x.len(),
            k_y.len()
        )));
    }
    for i in 0..n {
        let triplet = kernel(x[i], y, tier);
        (k[i], k_x[i], k_y[i]) = (triplet.k, triplet.k_x, triplet.k_y);
    }
    Ok(())
}
```

- [ ] **Step 5: Run the tests and lint checks**

Run:
```bash
cargo test -p absrbrs -- --nocapture
cargo clippy --workspace --all-targets -- -D warnings
cargo fmt --all -- --check
RUSTDOCFLAGS="-D warnings" cargo doc -p absrbrs --no-deps --document-private-items
```
Expected:
- every test passes;
- the Accurate tier's printed worst value error is around 1e-13 or below, and its derivative error is ≤ 1e-12;
- clippy, fmt and rustdoc clean.

If the value error exceeds 2e-13, the port differs from the source. Compare the region it occurs in, line by line, with the C code. Do not loosen the tolerance.

- [ ] **Step 6: Commit**

```bash
git add absrbrs/src absrbrs/tests/faddeeva_accuracy.rs
git commit -m "feat(faddeeva): Accurate tier w(z) with derivative; public API (w, kernel, kernel_into)

Co-Authored-By: <your model> <noreply@anthropic.com>"
```

---

### Task 3: Fast tier

**Files:**
- Create: `reference/absrbrs_reference/generators/weideman.py`, `reference/tests/test_weideman.py`
- Modify: `reference/absrbrs_reference/generators/__init__.py`, adding `"weideman": weideman.generate`
- Generate: `reference/data/weideman.jsonl` and the updated `reference/data/manifest.json`
- Create: `absrbrs/src/faddeeva/fast/mod.rs`, `absrbrs/src/faddeeva/fast/weideman.rs`, `absrbrs/src/faddeeva/fast/taylor.rs`
- Modify: `absrbrs/src/faddeeva/mod.rs` (`Tier::Fast`, default and dispatch), `absrbrs/Cargo.toml` (proptest dev-dependency, `exclude`)
- Test: `absrbrs/tests/faddeeva_accuracy.rs` (Fast cases), `absrbrs/tests/faddeeva_fast_vs_accurate.rs`

**Interfaces:**
- Consumes: from Task 2, `cf::laplace_cf`, `cf::gauss_cf`, `accurate::w_im`, `accurate::w_and_derivative`, `FRAC_1_SQRT_PI` and the public API.
- Produces:
  - `Tier::Fast`, which becomes `#[default]` (so add `Default` to the derive list);
  - `pub(crate) fn fast::w_and_derivative(z: Complex64) -> (Complex64, Complex64)` for Re z ≥ 0 and Im z ≥ 0;
  - `pub(crate) const fast::weideman::{N, L, A}`;
  - the reference file `weideman.jsonl`, with one record: inputs `{"n_terms": 24}`, outputs `L` (one value) and `a` (24 values, a₁ … a₂₄).

The Fast composite, from `docs/plans/phase-0-followups.md`. It works on x ≥ 0 and y ≥ 0, after the reflection in `w`. The first matching row wins, and r = |z| = `x.hypot(y)`:

| Order | Region | Method |
|---|---|---|
| 1 | y ≤ 0.01 and 6 ≤ x < 27 | `gauss_cf(z, 5)` |
| 2 | r ≥ 1750 | depth 0 (Lorentzian), in real arithmetic |
| 3 | 50 ≤ r < 1750 | depth 1, in real arithmetic |
| 4 | 15 ≤ r < 50 | `laplace_cf(z, 2)` |
| 5 | 9 ≤ r < 15 | `laplace_cf(z, 3)` |
| 6 | 8 ≤ r < 9 | `laplace_cf(z, 4)` |
| 7 | y ≤ 0.01 and x < 6 | Taylor series of order 3 |
| 8 | otherwise | Weideman, 24 terms |

- [ ] **Step 1: Generate the Weideman coefficients with Arb**

Create `reference/absrbrs_reference/generators/weideman.py`:
```python
"""Weideman's (1994) N = 24 rational approximation of w(z): the constant L and coefficients aₙ.

w(z) ≈ 2 Σₙ aₙ Zⁿ⁻¹/(L − iz)² + 1/(√π (L − iz)), Z = (L + iz)/(L − iz), L = √(N/√2).

The aₙ are the discrete Fourier coefficients of F(t) = exp(−t²)(L² + t²) at the 2M = 4N points
t_k = L tan(kπ/(2M)): aₙ = (F(0) + 2 Σₖ₌₁^(M−1) F(t_k) cos(πkn/M)) / (2M), n = 1 … N. (F(t_M) is
taken as 0, the limit as t → ∞.) This is the formula of the paper; its code listing is not used.
"""

from __future__ import annotations

from flint import arb, ctx

from absrbrs_reference.records import record

N = 24
PRECISION_BITS = 256


def coefficients() -> tuple[arb, list[arb]]:
    """L and [a₁, …, a_N] as Arb balls."""
    m = 2 * N
    big_l = (arb(N) / arb(2).sqrt()).sqrt()

    def f(k: int) -> arb:
        t = big_l * (arb.pi() * k / (2 * m)).tan()
        return (-(t * t)).exp() * (big_l * big_l + t * t)

    values = [f(k) for k in range(m)]
    a = [
        (values[0] + 2 * sum(values[k] * (arb.pi() * k * n / m).cos() for k in range(1, m)))
        / (2 * m)
        for n in range(1, N + 1)
    ]
    return big_l, a


def generate() -> tuple[int, list[dict]]:
    with ctx.workprec(PRECISION_BITS):
        big_l, a = coefficients()
        return PRECISION_BITS, [
            record("weideman", "n_terms=24", {"n_terms": N}, {"L": [big_l], "a": a})
        ]
```

Register it in `reference/absrbrs_reference/generators/__init__.py`, keeping the names sorted:
```python
from absrbrs_reference.generators import faddeeva, weideman
...
GENERATORS: dict[str, Generator] = {
    "faddeeva": faddeeva.generate,
    "weideman": weideman.generate,
}
```

Create `reference/tests/test_weideman.py`:
```python
"""The Weideman coefficients: certified, and copied exactly into the Rust Fast tier."""

import re
from pathlib import Path

from absrbrs_reference.generators import weideman

RUST = (
    Path(__file__).resolve().parents[2] / "absrbrs" / "src" / "faddeeva" / "fast" / "weideman.rs"
)


def test_values_are_tight_and_decay():
    _, records = weideman.generate()
    (r,) = records
    a = r["outputs"]["a"]
    assert len(a["hi"]) == 24
    assert all(e <= 1e-60 for e in a["err"])
    assert abs(a["hi"][23]) < 1e-9 < abs(a["hi"][0])


def test_rust_constants_are_the_nearest_doubles():
    _, (r,) = weideman.generate()
    text = RUST.read_text("utf-8")
    l_rust = float(re.search(r"pub\(crate\) const L: f64 = ([^;]+);", text).group(1).replace("_", ""))
    a_block = re.search(r"pub\(crate\) const A: \[f64; N\] = \[(.*?)\];", text, re.S).group(1)
    a_rust = [float(v.replace("_", "")) for v in re.findall(r"-?[0-9][0-9_.e+-]*", a_block)]
    assert l_rust == r["outputs"]["L"]["hi"][0]
    assert a_rust == list(reversed(r["outputs"]["a"]["hi"]))  # highest power first
```

Run:
```bash
uv run --locked --directory reference python -m absrbrs_reference generate --only weideman
uv run --locked --directory reference python -m pytest -q tests/test_weideman.py
```
Expected: the generator prints `weideman.jsonl: 1 records`. `test_values_are_tight_and_decay` passes, and `test_rust_constants_are_the_nearest_doubles` fails with `FileNotFoundError`, because `weideman.rs` does not exist yet.

- [ ] **Step 2: Write the failing Rust tests**

Append to `absrbrs/tests/faddeeva_accuracy.rs`:
```rust
#[test]
fn fast_tier_meets_its_targets() {
    let (v, d) = worst(Tier::Fast);
    eprintln!("Fast: value {:e} at ({}, {}), derivative {:e} at ({}, {})", v.0, v.1, v.2, d.0, d.1, d.2);
    assert!(v.0 <= 1e-6, "value error {:e} at x = {}, y = {}", v.0, v.1, v.2);
    assert!(d.0 <= 1e-5, "derivative error {:e} at x = {}, y = {}", d.0, d.1, d.2);
}

#[test]
fn fast_is_the_default_tier() {
    assert_eq!(Tier::default(), Tier::Fast);
}
```
In `negative_x_follows_the_symmetries_exactly`, change `for tier in [Tier::Accurate]` to `for tier in [Tier::Fast, Tier::Accurate]`.

Create `absrbrs/tests/faddeeva_fast_vs_accurate.rs`:
```rust
//! The Fast tier against the Accurate tier (itself within 2e-13 of Arb): on random points, and
//! on both sides of every region boundary of the Fast composite.

use absrbrs::faddeeva::{self, Complex64, Tier};
use proptest::prelude::*;

/// Value and derivative errors of the Fast tier at z, with the Accurate tier as reference.
fn errors(z: Complex64) -> (f64, f64) {
    let (w, dw) = faddeeva::w(z, Tier::Fast);
    let (wr, dwr) = faddeeva::w(z, Tier::Accurate);
    let ve = if wr.re < f64::MIN_POSITIVE {
        if w.re.is_finite() { 0.0 } else { f64::INFINITY }
    } else {
        (w.re - wr.re).abs() / wr.re
    };
    let de = (dw.re - dwr.re).abs().max((dw.im - dwr.im).abs()) / dwr.norm();
    (if ve.is_nan() { f64::INFINITY } else { ve }, if de.is_nan() { f64::INFINITY } else { de })
}

fn x_values() -> impl Strategy<Value = f64> {
    prop_oneof![0.0..30.0f64, 0.0..12.0f64, (-8.0..8.0f64).prop_map(|e| 10f64.powf(e))]
}

fn y_values() -> impl Strategy<Value = f64> {
    prop_oneof![
        Just(0.0),
        (-4.0..0.0f64).prop_map(|e| 10f64.powf(e)),
        (-30.0..6.0f64).prop_map(|e| 10f64.powf(e)),
    ]
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(20_000))]

    #[test]
    fn fast_tier_is_within_its_targets(x in x_values(), y in y_values(), negative in any::<bool>()) {
        let x = if negative { -x } else { x };
        let (ve, de) = errors(Complex64::new(x, y));
        prop_assert!(ve <= 1e-6, "value error {ve:e} at ({x}, {y})");
        prop_assert!(de <= 1e-5, "derivative error {de:e} at ({x}, {y})");
    }
}

/// The next double above a positive v.
fn up(v: f64) -> f64 {
    f64::from_bits(v.to_bits() + 1)
}

/// The next double below a positive v.
fn down(v: f64) -> f64 {
    f64::from_bits(v.to_bits() - 1)
}

/// Points on both sides of every boundary of the composite: the |z| circles at all angles, and
/// the x and y thresholds of the near-axis regions.
#[test]
fn fast_tier_is_within_its_targets_at_every_boundary() {
    let mut points = Vec::new();
    for r in [8.0, 9.0, 15.0, 50.0, 1750.0f64] {
        for j in 0..=200 {
            let theta = std::f64::consts::FRAC_PI_2 * f64::from(j) / 200.0;
            for s in [1.0 - 1e-15, 1.0, 1.0 + 1e-15] {
                points.push((r * s * theta.cos(), r * s * theta.sin()));
            }
        }
    }
    for x in [6.0, 27.0f64] {
        for y in [0.0, 1e-30, 1e-10, 1e-4, 0.01, up(0.01)] {
            points.extend([(down(x), y), (x, y)]);
        }
    }
    for i in 0..=600 {
        let x = f64::from(i) * 0.05;
        points.extend([(x, 0.01), (x, up(0.01))]);
    }
    let mut worst = (0.0f64, 0.0f64);
    for (x, y) in points {
        let (ve, de) = errors(Complex64::new(x, y.max(0.0)));
        assert!(ve <= 1e-6, "value error {ve:e} at ({x}, {y})");
        assert!(de <= 1e-5, "derivative error {de:e} at ({x}, {y})");
        worst = (worst.0.max(ve), worst.1.max(de));
    }
    eprintln!("boundaries: worst value {:e}, derivative {:e}", worst.0, worst.1);
}
```

(`f64::next_up` is newer than the MSRV, and clippy's `incompatible_msrv` lint applies to tests too; hence the bit-level helpers.)

Add the proptest dev-dependency to `absrbrs/Cargo.toml`, after running the usask-foundation:dependency-audit skill for proptest 1 (MIT OR Apache-2.0, dev only, also used by scatrs):
```toml
proptest = "1"
```
Change the `exclude` line to `exclude = ["/tests/", "*.proptest-regressions"]`.

Run: `cargo test -p absrbrs --test faddeeva_accuracy --test faddeeva_fast_vs_accurate`

Expected: compile errors ``no variant `Fast` `` and ``no function `default` ``.

- [ ] **Step 3: Write the Fast tier**

Create `absrbrs/src/faddeeva/fast/weideman.rs`. The constants are the nearest doubles to the Arb values from Step 1, and `reference/tests/test_weideman.py` checks them:
```rust
//! Weideman's (1994) N = 24 rational approximation of w(z), with the derivative of the formula.

use num_complex::Complex64;

use crate::faddeeva::FRAC_1_SQRT_PI;

/// Number of terms.
pub(crate) const N: usize = 24;
/// L = √(N/√2).
pub(crate) const L: f64 = 4.119534287814235;
/// a_N, …, a₁: highest power of Z first, for Horner's rule (reference/data/weideman.jsonl).
pub(crate) const A: [f64; N] = [
    -1.5137461285014412e-10,
    4.9048204900217325e-09,
    1.3310454420371002e-09,
    -3.008282348388502e-08,
    -1.9122258888735165e-08,
    1.8738343454069238e-07,
    2.5682641335791855e-07,
    -1.085647579569039e-06,
    -3.0388931842808516e-06,
    4.139461724333319e-06,
    3.0471066082991645e-05,
    2.4331415462260513e-05,
    -0.00020748431511424012,
    -0.0007816642995623505,
    -0.0004936426901285409,
    0.006215006362949199,
    0.033723366855315975,
    0.10838723484566723,
    0.26549639598807706,
    0.5361139535729115,
    0.9257087138588673,
    1.394819673379119,
    1.8562864992055401,
    2.1978589365315413,
];

const I: Complex64 = Complex64::new(0.0, 1.0);

/// w(z) ≈ 2 p(Z)/(L − iz)² + 1/(√π (L − iz)), Z = (L + iz)/(L − iz), p(Z) = Σ aₙ Zⁿ⁻¹, and its
/// derivative.
pub(crate) fn w_and_derivative(z: Complex64) -> (Complex64, Complex64) {
    let d = L - I * z;
    let inv = d.inv();
    let big_z = (L + I * z) * inv;
    let (mut p, mut dp) = (Complex64::new(A[0], 0.0), Complex64::new(0.0, 0.0));
    for &a in &A[1..] {
        dp = dp * big_z + p;
        p = p * big_z + a;
    }
    let inv2 = inv * inv;
    let w = 2.0 * p * inv2 + FRAC_1_SQRT_PI * inv;
    let dz = 2.0 * I * L * inv2; // dZ/dz
    let dw = 2.0 * dp * dz * inv2 + 4.0 * I * p * inv2 * inv + I * FRAC_1_SQRT_PI * inv2;
    (w, dw)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The coefficients agree with an f64 evaluation of the same discrete Fourier formula.
    #[test]
    fn coefficients_follow_from_the_formula() {
        let m = 2 * N;
        let f = |k: usize| {
            let t = L * (std::f64::consts::PI * k as f64 / (2 * m) as f64).tan();
            (-t * t).exp() * (L * L + t * t)
        };
        for n in 1..=N {
            let sum: f64 = (1..m)
                .map(|k| f(k) * (std::f64::consts::PI * (k * n) as f64 / m as f64).cos())
                .sum();
            let a = (f(0) + 2.0 * sum) / (2 * m) as f64;
            assert!((a - A[N - n]).abs() <= 1e-13, "a_{n}: {a} vs {}", A[N - n]);
        }
        assert!((L - (N as f64 / std::f64::consts::SQRT_2).sqrt()).abs() <= 1e-15 * L);
    }
}
```

Create `absrbrs/src/faddeeva/fast/taylor.rs`:
```rust
//! The small-y Taylor series of w about the real axis, with the formula's own derivatives.

use num_complex::Complex64;

use crate::faddeeva::{FRAC_1_SQRT_PI, accurate};

/// Order of the series.
const ORDER: usize = 3;

/// w(x + iy) ≈ Σₖ (iy)ᵏ w⁽ᵏ⁾(x)/k!, k = 0 … ORDER, with w(x) = exp(−x²) + i w_im(x) and
/// w⁽ᵏ⁺¹⁾ = −2x w⁽ᵏ⁾ − 2k w⁽ᵏ⁻¹⁾. The series is not analytic in z, so the derivative returned is
/// K_x − i K_y of the formula: K_x = Re Σₖ (iy)ᵏ w⁽ᵏ⁺¹⁾/k! and K_y = Re(i Σₖ₍ₖ<ORDER₎ (iy)ᵏ w⁽ᵏ⁺¹⁾/k!).
/// For |x| < 6 only: w′(x) comes from the identity −2xw + 2i/√π, which amplifies the error of
/// w_im by about 2x² (spec §6.1; accepted because w_im is accurate to about 1e-15).
pub(crate) fn w_and_derivative(x: f64, y: f64) -> (Complex64, Complex64) {
    let mut d = [Complex64::new(0.0, 0.0); ORDER + 2];
    d[0] = Complex64::new((-x * x).exp(), accurate::w_im(x));
    d[1] = -2.0 * x * d[0] + Complex64::new(0.0, 2.0 * FRAC_1_SQRT_PI);
    for k in 1..=ORDER {
        d[k + 1] = -2.0 * x * d[k] - 2.0 * k as f64 * d[k - 1];
    }
    let iy = Complex64::new(0.0, y);
    let zero = Complex64::new(0.0, 0.0);
    let (mut w, mut d_x, mut d_y) = (zero, zero, zero);
    let mut c = Complex64::new(1.0, 0.0); // (iy)ᵏ/k!
    for k in 0..=ORDER {
        w += c * d[k];
        d_x += c * d[k + 1];
        if k < ORDER {
            d_y += c * d[k + 1];
        }
        c = c * iy / (k + 1) as f64;
    }
    let k_x = d_x.re;
    let k_y = (Complex64::new(0.0, 1.0) * d_y).re;
    (w, Complex64::new(k_x, -k_y))
}
```

Create `absrbrs/src/faddeeva/fast/mod.rs`. Depths 0 and 1 hold about 98% of the points of atmospheric line windows, so they are written in real arithmetic, without complex division:
```rust
//! The Fast tier (spec §6.3): the composite chosen in phase 0 (docs/plans/phase-0-followups.md).

mod taylor;
mod weideman;

use num_complex::Complex64;

use super::FRAC_1_SQRT_PI;
use super::cf::{gauss_cf, laplace_cf};

/// w(z) and w′(z) for Re z ≥ 0 and Im z ≥ 0; the first matching region of the composite wins.
pub(crate) fn w_and_derivative(z: Complex64) -> (Complex64, Complex64) {
    let (x, y) = (z.re, z.im);
    if y <= 0.01 && (6.0..27.0).contains(&x) {
        return gauss_cf(z, 5);
    }
    let r = x.hypot(y);
    if r >= 1750.0 {
        return lorentzian(x, y);
    }
    if r >= 50.0 {
        return humlicek_1(x, y);
    }
    if r >= 15.0 {
        return laplace_cf(z, 2);
    }
    if r >= 9.0 {
        return laplace_cf(z, 3);
    }
    if r >= 8.0 {
        return laplace_cf(z, 4);
    }
    if y <= 0.01 {
        return taylor::w_and_derivative(x, y);
    }
    weideman::w_and_derivative(z)
}

/// Depth 0: w = i/(√π z) = (y + ix)/(√π r²), w′ = −i/(√π z²) = −(2xy + i(x − y)(x + y))/(√π r⁴).
fn lorentzian(x: f64, y: f64) -> (Complex64, Complex64) {
    let r2 = x * x + y * y;
    let s = FRAC_1_SQRT_PI / r2;
    let s2 = s / r2;
    (
        Complex64::new(y * s, x * s),
        Complex64::new(-2.0 * x * y * s2, -(x - y) * (x + y) * s2),
    )
}

/// Depth 1 (Humlíček's region I): w = (i/√π) z/D with D = z² − ½, w′ = (i/√π)(−z² − ½)/D².
/// With a = (x − y)(x + y) − ½ and b = 2xy: Re w = y(x² + y² + ½)/(√π |D|²) and
/// Im w = x(x² + y² − ½)/(√π |D|²), without cancellation.
fn humlicek_1(x: f64, y: f64) -> (Complex64, Complex64) {
    let r2 = x * x + y * y;
    let a = (x - y) * (x + y) - 0.5;
    let b = 2.0 * x * y;
    let d2 = a * a + b * b; // |D|²
    let s = FRAC_1_SQRT_PI / d2;
    let w = Complex64::new(y * (r2 + 0.5) * s, x * (r2 - 0.5) * s);
    // w′ = (i/√π) n conj(D)² / |D|⁴, n = −z² − ½ = (−a − 1, −b).
    let (nr, ni) = (-a - 1.0, -b);
    let (cr, ci) = ((a - b) * (a + b), -2.0 * a * b); // conj(D)²
    let (pr, pi) = (nr * cr - ni * ci, nr * ci + ni * cr);
    let t = s / d2;
    (w, Complex64::new(-pi * t, pr * t))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The real-arithmetic forms equal the complex continued fraction of the same depth.
    #[test]
    fn real_forms_match_the_continued_fraction() {
        for (x, y) in [(50.0, 0.0), (60.0, 1e-3), (200.0, 37.0), (1749.0, 5.0), (1750.0, 0.0), (1e5, 1e4), (3e3, 2e3)] {
            let z = Complex64::new(x, y);
            for (depth, (w, dw)) in [(0, lorentzian(x, y)), (1, humlicek_1(x, y))] {
                let (wc, dwc) = laplace_cf(z, depth);
                assert!((w.re - wc.re).abs() <= 1e-14 * wc.re.abs().max(f64::MIN_POSITIVE), "K, depth {depth} at ({x}, {y})");
                assert!((w - wc).norm() <= 1e-14 * wc.norm(), "w, depth {depth} at ({x}, {y})");
                assert!((dw - dwc).norm() <= 1e-14 * dwc.norm(), "w′, depth {depth} at ({x}, {y})");
            }
        }
    }
}
```
At y = 0, Re w is exactly 0 in both forms, so the `K` assertion compares 0 with 0. Keep `.max(f64::MIN_POSITIVE)` so that it does not divide by zero.

Update `absrbrs/src/faddeeva/mod.rs`:
- add `mod fast;`;
- add `Default` to `Tier`'s derive list;
- add the variant, before `Accurate`:
  ```rust
      /// Value error ≤ 1e-6 relative to Re w; derivative error ≤ 1e-5 relative to |w′|. The
      /// region-switched composite of docs/plans/phase-0-followups.md. The default.
      #[default]
      Fast,
  ```
- add the dispatch arm `Tier::Fast => fast::w_and_derivative(z),` to `w`.

- [ ] **Step 4: Run all the tests**

Run:
```bash
cargo test -p absrbrs -- --nocapture
uv run --locked --directory reference python -m pytest -q
cargo clippy --workspace --all-targets -- -D warnings
cargo fmt --all -- --check
```
Expected:
- Every Rust test passes. The printed Fast worst value error is ≤ about 5e-7, and the derivative error is ≤ about 2e-6; phase 0 measured 4.9e-7 and 1.8e-6.
- The reference tests: `84 passed`.

If proptest finds a failure, it writes `absrbrs/tests/faddeeva_fast_vs_accurate.proptest-regressions`. Investigate the failing point against the region table above; do not widen tolerances.

- [ ] **Step 5: Commit**

```bash
git add absrbrs reference
git commit -m "feat(faddeeva): Fast tier (phase 0 composite) and Weideman-24 coefficients from Arb

Dependency audit: proptest 1 (MIT OR Apache-2.0; dev only) ...

Co-Authored-By: <your model> <noreply@anthropic.com>"
```

---

### Task 4: Benchmarks

**Files:**
- Modify: `absrbrs/Cargo.toml` (criterion dev-dependency, `[[bench]]`)
- Create: `absrbrs/benches/faddeeva.rs`, `docs/plans/phase-1a-followups.md`
- Modify: `.github/workflows/ci.yml`, adding a bench smoke-test step to `lint`

**Interfaces:**
- Consumes: `absrbrs::faddeeva::{kernel, kernel_into, Tier}`.
- Produces: the bench `faddeeva` and the measured numbers in `docs/plans/phase-1a-followups.md`.

- [ ] **Step 1: Add criterion**

Run the usask-foundation:dependency-audit skill for criterion 0.7 (MIT OR Apache-2.0; dev only; used by scatrs with the same features). Then add to `absrbrs/Cargo.toml`:
```toml
[dev-dependencies]
criterion = { version = "0.7", default-features = false, features = ["cargo_bench_support"] }

[[bench]]
name = "faddeeva"
harness = false
```

- [ ] **Step 2: Write the benchmark**

Create `absrbrs/benches/faddeeva.rs`. The 8-pole J-pole baseline is sasktran2's (MIT, Copyright (c) 2024 University of Saskatchewan Atmospheric Research Group), copied here for comparison only.
```rust
//! Throughput of the Faddeeva kernels: one point in each Fast-tier region, and the grid points
//! of representative line windows (the phase 0 workload: O₂ A-band and CO₂ 1.6 µm lines at the
//! surface, near 30 km and near 65 km, 0.005 cm⁻¹ grid, 25 cm⁻¹ window), against sasktran2's
//! 8-pole J-pole approximation.

use std::hint::black_box;

use absrbrs::faddeeva::{Complex64, Tier, kernel, kernel_into};
use criterion::{Criterion, Throughput, criterion_group, criterion_main};

/// sasktran2's J-pole constants (rust/sasktran2-rs/src/math/errorfunctions/optimized.rs):
/// Copyright (c) 2024 University of Saskatchewan Atmospheric Research Group, MIT licence.
const JPOLE_B: [Complex64; 4] = [
    Complex64::new(0.00383968430671409, -0.0119854387180615),
    Complex64::new(-0.321597857664957, -0.218883985607935),
    Complex64::new(2.55515264319988, 0.613958600684469),
    Complex64::new(-2.73739446984183, 5.69007914897806),
];
const JPOLE_C: [Complex64; 4] = [
    Complex64::new(2.51506776338386, -1.60713668042405),
    Complex64::new(-1.68985621846204, -1.66471695485661),
    Complex64::new(0.981465428659098, -1.70017951305004),
    Complex64::new(-0.322078795578047, -1.71891780447016),
];
const FRAC_1_SQRT_PI: f64 = 0.564_189_583_547_756_3;

/// The J-pole kernel triplet (sasktran2 computes the value and derivative this way). The other
/// four poles are −conj of these, with conj coefficients.
fn jpole_kernel(x: f64, y: f64) -> (f64, f64, f64) {
    let z = Complex64::new(x, y);
    let (mut s, mut ds) = (Complex64::new(0.0, 0.0), Complex64::new(0.0, 0.0));
    for (b, c) in JPOLE_B.iter().zip(JPOLE_C.iter()) {
        for (b, c) in [(*b, *c), (b.conj(), -c.conj())] {
            let inv = (z - c).inv();
            s += b * inv;
            ds += b * inv * inv;
        }
    }
    let i = Complex64::new(0.0, 1.0);
    let w = -i * FRAC_1_SQRT_PI * s;
    let dw = i * FRAC_1_SQRT_PI * ds;
    (w.re, dw.re, -dw.im)
}

const K_B: f64 = 1.380649e-23;
const N_A: f64 = 6.02214076e23;
const C: f64 = 299_792_458.0;

/// (name, centre cm⁻¹, mass g/mol, γ_air cm⁻¹/atm, n_air).
const LINES: [(&str, f64, f64, f64, f64); 2] =
    [("o2", 13100.0, 31.98983, 0.045, 0.72), ("co2", 6300.0, 43.98983, 0.070, 0.75)];
/// (name, temperature K, pressure Pa).
const CONDITIONS: [(&str, f64, f64); 3] =
    [("surface", 250.0, 101325.0), ("30km", 225.0, 1000.0), ("65km", 200.0, 10.0)];

/// x at every grid point of a line's window, and the line's y.
fn window(centre: f64, mass: f64, gamma: f64, n: f64, t: f64, p: f64) -> (Vec<f64>, f64) {
    let sigma_d = centre / C * (2.0 * K_B * t * N_A * 1e3 / mass).sqrt();
    let y = gamma * (296.0 / t).powf(n) * p / 101325.0 / sigma_d;
    ((-5000..=5000).map(|i| f64::from(i) * 0.005 / sigma_d).collect(), y)
}

fn bench_regions(c: &mut Criterion) {
    let mut group = c.benchmark_group("point");
    let points = [
        ("weideman", 2.0, 0.5),
        ("taylor", 3.0, 1e-3),
        ("exp_cf5", 10.0, 1e-3),
        ("cf4", 8.5, 0.5),
        ("cf3", 12.0, 1.0),
        ("cf2", 30.0, 1.0),
        ("cf1", 300.0, 1.0),
        ("cf0", 3000.0, 1.0),
    ];
    for (name, x, y) in points {
        group.bench_function(format!("fast/{name}"), |b| {
            b.iter(|| kernel(black_box(x), black_box(y), Tier::Fast))
        });
        group.bench_function(format!("accurate/{name}"), |b| {
            b.iter(|| kernel(black_box(x), black_box(y), Tier::Accurate))
        });
        group.bench_function(format!("jpole/{name}"), |b| {
            b.iter(|| jpole_kernel(black_box(x), black_box(y)))
        });
    }
    group.finish();
}

fn bench_lines(c: &mut Criterion) {
    let mut group = c.benchmark_group("line");
    for (line, centre, mass, gamma, n) in LINES {
        for (cond, t, p) in CONDITIONS {
            let (x, y) = window(centre, mass, gamma, n, t, p);
            let len = x.len();
            let (mut k, mut kx, mut ky) = (vec![0.0; len], vec![0.0; len], vec![0.0; len]);
            group.throughput(Throughput::Elements(len as u64));
            for (tier_name, tier) in [("fast", Tier::Fast), ("accurate", Tier::Accurate)] {
                group.bench_function(format!("{tier_name}/{line}_{cond}"), |b| {
                    b.iter(|| kernel_into(black_box(&x), y, tier, &mut k, &mut kx, &mut ky).unwrap())
                });
            }
            group.bench_function(format!("jpole/{line}_{cond}"), |b| {
                b.iter(|| {
                    for (i, &xi) in x.iter().enumerate() {
                        (k[i], kx[i], ky[i]) = jpole_kernel(black_box(xi), y);
                    }
                })
            });
        }
    }
    group.finish();
}

criterion_group!(benches, bench_regions, bench_lines);
criterion_main!(benches);
```

- [ ] **Step 3: Run it in test mode, then for real**

Run:
```bash
cargo test -p absrbrs --bench faddeeva --locked
cargo bench -p absrbrs --bench faddeeva
```
Expected:
- The test mode reports each benchmark as `Success`.
- `cargo bench` prints a time per point for every `point/…` benchmark, and a time and throughput (elements/s) for every `line/…` benchmark. It writes them to `target/criterion/`.

- [ ] **Step 4: Record the numbers and add the CI smoke test**

Create `docs/plans/phase-1a-followups.md`:
```markdown
# Phase 1a follow-ups and decisions

## Benchmarks (`cargo bench -p absrbrs --bench faddeeva`)

Measured on <machine: CPU model, OS>, <date>, release profile (fat LTO, one codegen unit).
Times are criterion's median estimates.

| Region | Fast (ns/point) | Accurate (ns/point) | J-pole (ns/point) |
|---|---|---|---|
| weideman | | | |
| taylor | | | |
| exp_cf5 | | | |
| cf4 | | | |
| cf3 | | | |
| cf2 | | | |
| cf1 | | | |
| cf0 | | | |

| Line window | Fast (ns/point) | Accurate (ns/point) | J-pole (ns/point) |
|---|---|---|---|
| o2_surface | | | |
| o2_30km | | | |
| o2_65km | | | |
| co2_surface | | | |
| co2_30km | | | |
| co2_65km | | | |

<Two or three sentences: the Fast tier's speed relative to the J-pole on line windows (the phase
0 cost proxy predicted about 4× fewer flop-equivalents), and whether the depth-0/1
real-arithmetic forms dominate as expected.>
```
Fill in every cell from the `cargo bench` output: for a line benchmark, the median time divided by 10 001 points. Fill in the sentences. Do not leave placeholders.

In `.github/workflows/ci.yml`'s `lint` job, after the step `Build the libraries`, add:
```yaml
      - name: Run each benchmark once (criterion's test mode; timing is not checked)
        run: cargo test -p absrbrs --bench faddeeva --locked
```
Add the same command to the "Checks" block of `CONTRIBUTING.md`, after `cargo test --workspace`. Also add the line ``Wall-clock benchmarks of the Faddeeva kernels: `cargo bench -p absrbrs --bench faddeeva`.`` below the block.

Run: `uvx --from actionlint-py==1.7.12.25 actionlint`. Expected: no output.

- [ ] **Step 5: Commit**

```bash
git add absrbrs/Cargo.toml Cargo.lock absrbrs/benches docs/plans/phase-1a-followups.md .github/workflows/ci.yml CONTRIBUTING.md
git commit -m "bench(faddeeva): kernels per region and over line windows, against sasktran2's J-pole

Dependency audit: criterion 0.7 (MIT OR Apache-2.0; dev only) ...

Co-Authored-By: <your model> <noreply@anthropic.com>"
```

---

### Task 5: Python bindings (`absrbrs.special`)

**Files:**
- Modify: `absrbrs-py/Cargo.toml` (numpy), `absrbrs-py/src/lib.rs`
- Create: `absrbrs-py/src/special.rs`, `python/absrbrs/special.py`
- Modify: `python/absrbrs/__init__.py`, `python/absrbrs/_absrbrs/__init__.pyi`, `tests/test_package.py`, `docs/api.md`
- Test: `tests/test_special.py`

**Interfaces:**
- Consumes: `absrbrs::faddeeva::{w, kernel, Tier, Complex64}` and `crate::errors::to_pyerr`.
- Produces:
  - in the extension: `faddeeva_w(z: NDArray[complex128] (1-D, contiguous), tier: str) -> tuple[NDArray[complex128], NDArray[complex128]]` and `voigt_kernel(x: NDArray[float64], y: NDArray[float64], tier: str) -> tuple[NDArray[float64], NDArray[float64], NDArray[float64]]`;
  - in Python: `absrbrs.special.faddeeva(z, *, tier="fast") -> (w, dw)` and `absrbrs.special.voigt_kernel(x, y, *, tier="fast") -> (k, k_x, k_y)`, both broadcasting and preserving shape.

- [ ] **Step 1: Write the failing tests**

Create `tests/test_special.py`:
```python
import numpy as np
import pytest

import absrbrs
from absrbrs import special

SQRT_PI = np.sqrt(np.pi)


def test_known_values():
    for tier in ("fast", "accurate"):
        w, dw = special.faddeeva([0.0, 1.0j, 1.0], tier=tier)
        np.testing.assert_allclose(w.real, [1.0, 0.42758357615580700, np.exp(-1.0)], rtol=1e-6)
        np.testing.assert_allclose(w.imag, [0.0, 0.0, 0.60715770584139372], rtol=1e-6, atol=1e-15)
        np.testing.assert_allclose(dw[0], 2j / SQRT_PI, rtol=1e-6)


def test_tiers_agree_within_the_fast_targets():
    rng = np.random.default_rng(1)
    z = rng.uniform(0, 30, 2000) + 1j * 10.0 ** rng.uniform(-6, 2, 2000)
    w_f, dw_f = special.faddeeva(z, tier="fast")
    w_a, dw_a = special.faddeeva(z, tier="accurate")
    assert np.max(np.abs(w_f.real - w_a.real) / w_a.real) <= 1e-6
    assert np.max(np.abs(dw_f - dw_a) / np.abs(dw_a)) <= 1e-5


def test_fast_is_the_default():
    z = np.array([0.3 + 0.2j, 12.0 + 1e-3j])
    for got, expected in zip(special.faddeeva(z), special.faddeeva(z, tier="fast"), strict=True):
        np.testing.assert_array_equal(got, expected)


def test_shapes_broadcast_and_scalars_give_zero_d_arrays():
    w, dw = special.faddeeva(np.ones((2, 3)) * (1 + 1j))
    assert w.shape == dw.shape == (2, 3)
    w, _ = special.faddeeva(1.0)
    assert w.shape == ()
    k, k_x, k_y = special.voigt_kernel(np.linspace(0, 5, 4)[:, None], [0.0, 0.5, 2.0])
    assert k.shape == k_x.shape == k_y.shape == (4, 3)


def test_voigt_kernel_is_the_triplet_of_w():
    x, y = np.linspace(-20, 20, 81), 0.25
    k, k_x, k_y = special.voigt_kernel(x, y)
    w, dw = special.faddeeva(x + 1j * y)
    np.testing.assert_array_equal(k, w.real)
    np.testing.assert_array_equal(k_x, dw.real)
    np.testing.assert_array_equal(k_y, -dw.imag)


def test_profile_normalization():
    # ∫ K dx = √π for any y.
    x = np.linspace(-2000, 2000, 400_001)
    k, _, _ = special.voigt_kernel(x, 1.0)
    integral = np.sum((k[1:] + k[:-1]) * np.diff(x)) / 2  # trapezoid rule (numpy 1.26 has no trapezoid)
    assert integral == pytest.approx(SQRT_PI, rel=1e-3)


@pytest.mark.parametrize(
    "z", [1.0 - 1e-3j, complex(1.0, np.nan), complex(np.inf, 1.0)], ids=["negative", "nan", "inf"]
)
def test_invalid_points_raise(z):
    with pytest.raises(absrbrs.InvalidInputError, match="Im z"):
        special.faddeeva([0.0, z])


def test_invalid_kernel_points_raise():
    with pytest.raises(absrbrs.InvalidInputError):
        special.voigt_kernel(1.0, -1.0)


def test_unknown_tier_raises():
    with pytest.raises(absrbrs.InvalidInputError, match="tier"):
        special.faddeeva(1.0, tier="exact")


def test_non_numeric_input_raises_type_error():
    with pytest.raises(TypeError):
        special.faddeeva("1.0")
```

In `tests/test_package.py`, update `test_test_hooks_live_outside_the_extension_namespace`'s expected `_absrbrs.__all__` to `{"__version__", "faddeeva_w", "voigt_kernel"}`. (`test_public_names` iterates `absrbrs.__all__`, so it covers `special` once Step 3 adds it there.)

Run: `uv run --locked pytest -q tests/test_special.py`

Expected: `ImportError: cannot import name 'special' from 'absrbrs'`.

- [ ] **Step 2: Write the extension functions**

Run the usask-foundation:dependency-audit skill for the numpy crate 0.29 (BSD-2-Clause; the bindings only; used by scatrs). Add to `absrbrs-py/Cargo.toml` `[dependencies]`:
```toml
numpy = "0.29"
```

Create `absrbrs-py/src/special.rs`:
```rust
//! Bindings of `absrbrs::faddeeva`, for `absrbrs.special`. Inputs are 1-D contiguous arrays
//! (the Python layer broadcasts and flattens); values are validated here.

use absrbrs::faddeeva::{self, Complex64, Tier};
use numpy::{IntoPyArray, PyArray1, PyArrayMethods};
use pyo3::exceptions::PyBufferError;
use pyo3::prelude::*;

use crate::errors::to_pyerr;

fn tier(py: Python<'_>, name: &str) -> PyResult<Tier> {
    match name {
        "fast" => Ok(Tier::Fast),
        "accurate" => Ok(Tier::Accurate),
        _ => Err(to_pyerr(
            py,
            absrbrs::Error::InvalidInput(format!("tier must be 'fast' or 'accurate', not {name:?}")),
        )),
    }
}

fn invalid(py: Python<'_>, message: String) -> PyErr {
    to_pyerr(py, absrbrs::Error::InvalidInput(message))
}

/// w(z) and w′(z) at every z (Im z ≥ 0, finite).
#[pyfunction]
pub fn faddeeva_w<'py>(
    py: Python<'py>,
    z: Bound<'py, PyArray1<Complex64>>,
    tier: &str,
) -> PyResult<(Bound<'py, PyArray1<Complex64>>, Bound<'py, PyArray1<Complex64>>)> {
    let tier = self::tier(py, tier)?;
    let z = z.try_readonly().map_err(|e| PyBufferError::new_err(e.to_string()))?;
    let z = z.as_slice()?;
    if let Some(bad) = z.iter().find(|z| !(z.im >= 0.0 && z.re.is_finite() && z.im.is_finite())) {
        return Err(invalid(py, format!("Im z must be >= 0 and z finite, got {bad}")));
    }
    let (w, dw): (Vec<_>, Vec<_>) = py.detach(|| z.iter().map(|&z| faddeeva::w(z, tier)).unzip());
    Ok((w.into_pyarray(py), dw.into_pyarray(py)))
}

/// The kernel triplet (K, K_x, K_y) at every (x, y) pair (y ≥ 0, both finite).
#[pyfunction]
#[allow(clippy::type_complexity)]
pub fn voigt_kernel<'py>(
    py: Python<'py>,
    x: Bound<'py, PyArray1<f64>>,
    y: Bound<'py, PyArray1<f64>>,
    tier: &str,
) -> PyResult<(Bound<'py, PyArray1<f64>>, Bound<'py, PyArray1<f64>>, Bound<'py, PyArray1<f64>>)> {
    let tier = self::tier(py, tier)?;
    let x = x.try_readonly().map_err(|e| PyBufferError::new_err(e.to_string()))?;
    let y = y.try_readonly().map_err(|e| PyBufferError::new_err(e.to_string()))?;
    let (x, y) = (x.as_slice()?, y.as_slice()?);
    if x.len() != y.len() {
        return Err(invalid(py, format!("x and y have {} and {} values", x.len(), y.len())));
    }
    if let Some((a, b)) = x.iter().zip(y).find(|(a, b)| !(a.is_finite() && b.is_finite() && **b >= 0.0)) {
        return Err(invalid(py, format!("y must be >= 0 and x, y finite, got x = {a}, y = {b}")));
    }
    let (k, (k_x, k_y)): (Vec<_>, (Vec<_>, Vec<_>)) = py.detach(|| {
        x.iter()
            .zip(y)
            .map(|(&x, &y)| {
                let t = faddeeva::kernel(x, y, tier);
                (t.k, (t.k_x, t.k_y))
            })
            .unzip()
    });
    Ok((k.into_pyarray(py), k_x.into_pyarray(py), k_y.into_pyarray(py)))
}
```

In `absrbrs-py/src/lib.rs`, add `mod special;` and, in `_absrbrs`, before `testing::register(m)?;`:
```rust
    m.add_function(wrap_pyfunction!(special::faddeeva_w, m)?)?;
    m.add_function(wrap_pyfunction!(special::voigt_kernel, m)?)?;
```

- [ ] **Step 3: Write the Python module and stubs**

Create `python/absrbrs/special.py`:
```python
"""Special functions: the Faddeeva function and the Voigt kernel (spec §6).

Both accuracy tiers are available: ``"fast"`` (the default; value error ≤ 1e-6 relative to
Re w, derivative error ≤ 1e-5 relative to |w′|) and ``"accurate"`` (≤ 2e-13 and ≤ 2e-12).
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from absrbrs import _absrbrs

__all__ = ["faddeeva", "voigt_kernel"]


def _numeric(name: str, value: ArrayLike) -> NDArray:
    array = np.asarray(value)
    if array.dtype.kind not in "biufc":
        raise TypeError(f"{name} must be numeric, not {array.dtype}")
    return array


def faddeeva(
    z: ArrayLike, *, tier: str = "fast"
) -> tuple[NDArray[np.complex128], NDArray[np.complex128]]:
    """The Faddeeva function w(z) = exp(−z²) erfc(−iz) and its derivative w′(z).

    Parameters
    ----------
    z : array_like
        Points with Im z ≥ 0 (complex or real).
    tier : {"fast", "accurate"}, optional
        The accuracy tier.

    Returns
    -------
    w, dw : numpy.ndarray
        w(z) and w′(z), complex, shaped like ``z``.

    Raises
    ------
    InvalidInputError
        If any Im z < 0, any z is not finite, or ``tier`` is unknown.
    """
    array = np.ascontiguousarray(_numeric("z", z), dtype=np.complex128)
    w, dw = _absrbrs.faddeeva_w(array.ravel(), tier)
    return w.reshape(array.shape), dw.reshape(array.shape)


def voigt_kernel(
    x: ArrayLike, y: ArrayLike, *, tier: str = "fast"
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    """The Voigt kernel triplet K = Re w(x + iy), K_x = ∂K/∂x and K_y = ∂K/∂y.

    Parameters
    ----------
    x, y : array_like
        Real points, broadcast together; y ≥ 0.
    tier : {"fast", "accurate"}, optional
        The accuracy tier.

    Returns
    -------
    k, k_x, k_y : numpy.ndarray
        Shaped like the broadcast of ``x`` and ``y``.

    Raises
    ------
    InvalidInputError
        If any y < 0, any value is not finite, or ``tier`` is unknown.
    """
    xb, yb = np.broadcast_arrays(
        _numeric("x", x).astype(np.float64), _numeric("y", y).astype(np.float64)
    )
    shape = xb.shape
    k, k_x, k_y = _absrbrs.voigt_kernel(
        np.ascontiguousarray(xb).ravel(), np.ascontiguousarray(yb).ravel(), tier
    )
    return k.reshape(shape), k_x.reshape(shape), k_y.reshape(shape)
```

In `python/absrbrs/__init__.py`, change `from absrbrs import _absrbrs` to `from absrbrs import _absrbrs, special`, and add `"special"` to `__all__` (sorted).

In `python/absrbrs/_absrbrs/__init__.pyi`, add `import numpy as np` and `from numpy.typing import NDArray`, extend `__all__` to `["__version__", "faddeeva_w", "voigt_kernel"]`, and add:
```python
def faddeeva_w(
    z: NDArray[np.complex128], tier: str
) -> tuple[NDArray[np.complex128], NDArray[np.complex128]]: ...
def voigt_kernel(
    x: NDArray[np.float64], y: NDArray[np.float64], tier: str
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]: ...
```

In `docs/api.md`, add before `## \`absrbrs.errors\``:
````markdown
## `absrbrs.special`

```{eval-rst}
.. automodule:: absrbrs.special
   :members:
```
````
Change the sentence "The package also re-exports the exceptions of {mod}`absrbrs.errors`." to "The package also re-exports the exceptions of {mod}`absrbrs.errors`, and provides {mod}`absrbrs.special`."

- [ ] **Step 4: Run every Python check**

Run:
```bash
uv sync
uv run --locked pytest -q
uv run --locked python -m mypy.stubtest absrbrs._absrbrs
uv run --locked ruff check . && uv run --locked ruff format --check .
cargo clippy --workspace --all-targets -- -D warnings
uv run --locked --extra docs sphinx-build -E -W --keep-going -b html docs docs/_build/html
```
Expected:
- pytest: `32 passed` (20 earlier and 12 new);
- stubtest, ruff and clippy clean;
- the Sphinx build succeeds with no warnings.

If numpydoc reports an unresolved type, add an alias to `numpydoc_xref_aliases` in `docs/conf.py`, e.g. `"numpy.ndarray": ":class:`numpy.ndarray`"`.

- [ ] **Step 5: Commit**

```bash
git add absrbrs-py python tests docs/api.md docs/conf.py Cargo.lock uv.lock
git commit -m "feat(python): absrbrs.special.faddeeva and voigt_kernel, both tiers

Dependency audit: numpy (Rust crate) 0.29 (BSD-2-Clause) ...

Co-Authored-By: <your model> <noreply@anthropic.com>"
```

---

### Task 6: Validation against Arb at CI time, documentation and wrap-up

**Files:**
- Create: `validation/pyproject.toml`, `validation/.python-version`, `validation/uv.lock` (generated)
- Create: `validation/absrbrs_validation/{__init__.py, arb.py, faddeeva.py, __main__.py}`, `validation/tests/test_faddeeva.py`
- Modify: `.gitignore`, `.github/workflows/ci.yml` (new `validation` job, in `all-checks`), `.github/dependabot.yml`, `pyproject.toml` (ruff `src`), `CONTRIBUTING.md`
- Generate: `docs/validation/data/faddeeva-accuracy.json`, `docs/validation/faddeeva-accuracy-table.md`
- Create: `docs/theory/faddeeva.md`, `docs/validation/faddeeva-accuracy.md`
- Modify: `docs/index.md`, `docs/plans/phase-1a-followups.md`, `CHANGELOG.md`

**Interfaces:**
- Consumes: `absrbrs.special.faddeeva`, through the installed package.
- Produces: the CI job `validation` and the generated accuracy report.

- [ ] **Step 1: Create the project**

Create `validation/pyproject.toml`:
```toml
[project]
name = "absrbrs-validation"
version = "0.1.0"
description = "Accuracy of absrbrs against Arb references computed at run time"
requires-python = ">=3.11"
# python-flint is pinned exactly, as in reference/.
dependencies = ["absrbrs", "numpy>=2.0", "python-flint==0.9.0"]

[dependency-groups]
dev = ["pytest>=8"]

[tool.uv]
package = false

[tool.uv.sources]
absrbrs = { path = "..", editable = true }

[tool.pytest.ini_options]
pythonpath = ["."]
testpaths = ["tests"]
```
```bash
mkdir -p validation/absrbrs_validation validation/tests
echo 3.12 > validation/.python-version
echo '"""Accuracy validation of absrbrs against Arb, computed at run time."""' > validation/absrbrs_validation/__init__.py
echo "/validation/.cache/" >> .gitignore
```
In the root `pyproject.toml`, add `"validation"` to `[tool.ruff] src`, after `"prototypes/faddeeva"`.

- [ ] **Step 2: Write the Arb evaluator and the measurement**

Create `validation/absrbrs_validation/arb.py`:
```python
"""Arb values of w(z) and w′(z), in parallel (as in prototypes/faddeeva, without its cache)."""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor

import numpy as np
from flint import acb, arb, ctx
from numpy.typing import NDArray

REL_BITS = 60


def arb_w(xy: tuple[float, float]) -> tuple[float, float, float, float]:
    """(Re w, Im w, Re w′, Im w′) at x + iy.

    Re w is known to 2^-60 relative, and each part of w′ to 2^-60 |w′|.
    """
    x, y = xy
    eps = arb(2) ** -REL_BITS
    bits = 128
    while bits <= 16384:
        with ctx.workprec(bits):
            z = acb(x, y)
            w = (-(z * z)).exp() * (-acb(0, 1) * z).erfc()
            dw = -2 * z * w + 2 * acb(0, 1) / arb.pi().sqrt()
            scale = abs(dw).mid()
            if (
                w.real.rad() <= abs(w.real.mid()) * eps
                and dw.real.rad() <= scale * eps
                and dw.imag.rad() <= scale * eps
            ):
                return (
                    float(w.real.mid()),
                    float(w.imag.mid()),
                    float(dw.real.mid()),
                    float(dw.imag.mid()),
                )
        bits *= 2
    raise ArithmeticError(f"w({x!r} + {y!r}i) did not reach 2^-{REL_BITS} by 16384 bits")


def compute(z: NDArray[np.complex128]) -> tuple[NDArray[np.complex128], NDArray[np.complex128]]:
    """Arb values (w, dw) at every point of z."""
    pairs = list(zip(z.real.tolist(), z.imag.tolist(), strict=True))
    with ProcessPoolExecutor() as pool:
        out = np.array(list(pool.map(arb_w, pairs, chunksize=512)))
    return out[:, 0] + 1j * out[:, 1], out[:, 2] + 1j * out[:, 3]
```

Create `validation/absrbrs_validation/faddeeva.py`:
```python
"""Both Faddeeva tiers against Arb: a dense grid (with every Fast-tier boundary and the double
on its other side) and random points, computed at run time (spec §6.2 metrics)."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from absrbrs import special
from absrbrs_validation import arb

TARGETS = {"fast": (1e-6, 1e-5), "accurate": (2e-13, 2e-12)}
X_EDGES = np.array([6.0, 8.0, 9.0, 15.0, 27.0, 50.0, 1750.0])
Y_EDGES = np.array([0.01, 0.1])
RANDOM_POINTS = 200_000
SEED = 20261015


def points() -> NDArray[np.complex128]:
    """The dense grid and the random points, x ≥ 0 and y ≥ 0."""
    xs = np.unique(
        np.concatenate(
            [
                np.arange(385) / 32.0,
                10.0 ** (np.arange(-128, 129) / 16),
                X_EDGES,
                np.nextafter(X_EDGES, -np.inf),
            ]
        )
    )
    ys = np.unique(
        np.concatenate(
            [[0.0], 10.0 ** (np.arange(-180, 37) / 6), Y_EDGES, np.nextafter(Y_EDGES, np.inf)]
        )
    )
    gx, gy = np.meshgrid(xs, ys)
    rng = np.random.default_rng(SEED)
    rx = np.where(
        rng.integers(0, 2, RANDOM_POINTS) == 0,
        rng.uniform(0, 30, RANDOM_POINTS),
        10.0 ** rng.uniform(-8, 8, RANDOM_POINTS),
    )
    ry = np.where(
        rng.integers(0, 10, RANDOM_POINTS) == 0, 0.0, 10.0 ** rng.uniform(-30, 6, RANDOM_POINTS)
    )
    x = np.concatenate([gx.ravel(), rx])
    y = np.concatenate([gy.ravel(), ry])
    z = np.empty(x.shape, dtype=np.complex128)
    z.real, z.imag = x, y
    return z


def measure(z: NDArray[np.complex128], w_ref, dw_ref, tier: str) -> dict:
    """The largest value and derivative errors of ``tier``, where they occur, and how many points
    were exempt (true K below the smallest normal double)."""
    w, dw = special.faddeeva(z, tier=tier)
    k = w_ref.real
    tiny = np.finfo(np.float64).tiny
    with np.errstate(all="ignore"):
        ve = np.where(k < tiny, np.where(np.isfinite(w.real), 0.0, np.inf), np.abs(w.real - k) / k)
        de = np.maximum(np.abs(dw.real - dw_ref.real), np.abs(dw.imag - dw_ref.imag)) / np.abs(dw_ref)
    ve = np.where(np.isnan(ve), np.inf, ve)
    de = np.where(np.isnan(de), np.inf, de)
    iv, idv = int(np.argmax(ve)), int(np.argmax(de))
    return {
        "tier": tier,
        "points": int(z.size),
        "exempt": int(np.count_nonzero(k < tiny)),
        "max_value_error": float(ve[iv]),
        "at_value": [float(z[iv].real), float(z[iv].imag)],
        "max_derivative_error": float(de[idv]),
        "at_derivative": [float(z[idv].real), float(z[idv].imag)],
        "value_target": TARGETS[tier][0],
        "derivative_target": TARGETS[tier][1],
    }
```

Create `validation/absrbrs_validation/__main__.py`:
```python
"""``python -m absrbrs_validation faddeeva``: measure both tiers, write the report and table."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from absrbrs_validation import arb, faddeeva

DOCS = Path(__file__).resolve().parents[2] / "docs" / "validation"


def table(results: list[dict]) -> str:
    lines = [
        "| Tier | Points | Max value error | Target | Max derivative error | Target |",
        "|---|---|---|---|---|---|",
    ]
    for r in results:
        lines.append(
            f"| {r['tier']} | {r['points']} | {r['max_value_error']:.1e} | {r['value_target']:.0e} "
            f"| {r['max_derivative_error']:.1e} | {r['derivative_target']:.0e} |"
        )
    return "\n".join(lines) + "\n"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="absrbrs_validation")
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("faddeeva", help="both Faddeeva tiers against Arb")
    parser.parse_args(argv)
    z = faddeeva.points()
    w, dw = arb.compute(z)
    results = [faddeeva.measure(z, w, dw, tier) for tier in ("fast", "accurate")]
    (DOCS / "data").mkdir(parents=True, exist_ok=True)
    (DOCS / "data" / "faddeeva-accuracy.json").write_text(
        json.dumps(results, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8"
    )
    (DOCS / "faddeeva-accuracy-table.md").write_text(table(results), encoding="utf-8")
    print(table(results), file=sys.stderr)
    ok = all(
        r["max_value_error"] <= r["value_target"]
        and r["max_derivative_error"] <= r["derivative_target"]
        for r in results
    )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
```

- [ ] **Step 3: Write the tests**

Create `validation/tests/test_faddeeva.py`:
```python
"""Both tiers meet their spec §6.3 targets on about 340 000 Arb points computed now."""

import pytest

from absrbrs_validation import arb, faddeeva


@pytest.fixture(scope="session")
def reference():
    z = faddeeva.points()
    w, dw = arb.compute(z)
    return z, w, dw


@pytest.mark.parametrize("tier", ["fast", "accurate"])
def test_tier_meets_its_targets(reference, tier):
    r = faddeeva.measure(*reference, tier)
    assert r["max_value_error"] <= r["value_target"], r
    assert r["max_derivative_error"] <= r["derivative_target"], r


def test_points_cover_the_domain_and_the_boundaries():
    z = faddeeva.points()
    assert z.real.min() == 0.0 and z.real.max() == 1e8
    assert z.imag.min() == 0.0 and z.imag.max() == pytest.approx(1e6)
    for edge in faddeeva.X_EDGES:
        assert (z.real == edge).any()
```

Run:
```bash
uv lock --directory validation
uv run --locked --directory validation python -m pytest -q
uv run --locked --directory validation python -m absrbrs_validation faddeeva
```
Expected:
- pytest: `3 passed`. Computing Arb at about 340 000 points takes a few seconds with all cores.
- The command exits with status 0, prints the table and writes the two docs files.
- The Fast tier's errors are around 5e-7 (value) and 2e-6 (derivative). The Accurate tier's are around 1e-13 and 1e-12.

- [ ] **Step 4: Add it to CI and dependabot**

In `.github/workflows/ci.yml`, add a job after `reference`:
```yaml
  validation:
    name: Accuracy against Arb (computed at run time)
    runs-on: ubuntu-latest
    timeout-minutes: 30
    env:
      CARGO_PROFILE_RELEASE_LTO: "false"
      CARGO_PROFILE_RELEASE_CODEGEN_UNITS: "16"
    steps:
      - uses: actions/checkout@3d3c42e5aac5ba805825da76410c181273ba90b1 # v7.0.1
        with:
          persist-credentials: false
      - name: Install Rust
        run: rustup toolchain install stable --profile minimal --no-self-update
      - uses: Swatinem/rust-cache@6323deb102c322ba6fcbdcafc7e3dddab59af2b6 # v2.9.2
        with:
          save-if: ${{ github.ref == 'refs/heads/main' }}
      - uses: astral-sh/setup-uv@1c37ad07a6a961277cf70c0d37d6f313000f5884 # v10.3.0
      - name: Both Faddeeva tiers against Arb (dense grid and random points)
        run: uv run --locked --directory validation python -m pytest -q
```
Add `validation` to `all-checks.needs`: `[lint, rust-test, msrv, python, docs, reference, validation]`.

In `.github/dependabot.yml`, add an entry for `/validation`, mirroring `/reference`, with `ignore: [{dependency-name: python-flint}]` and the comment "python-flint is pinned exactly".

In `CONTRIBUTING.md`'s "Checks" block, add `uv run --directory validation python -m pytest`. Under "Continuous integration", add the bullet "the accuracy of both Faddeeva tiers against Arb values computed at run time (`validation/`)".

Run: `uvx --from actionlint-py==1.7.12.25 actionlint`. Expected: no output.

- [ ] **Step 5: Write the documentation pages**

Create `docs/theory/faddeeva.md`:
```markdown
# The Faddeeva function and the Voigt kernel

Every Voigt line shape is K(x, y) = Re w(x + iy), the real part of the Faddeeva function
w(z) = exp(−z²) erfc(−iz), at x = (ν − ν₀′)/σ_D and y = γ_L/σ_D ≥ 0 (spec §5–6). The Jacobians
need its partial derivatives: K_x = Re w′ and K_y = −Im w′, since w is analytic. absrbrs computes
the triplet (K, K_x, K_y) at two accuracy tiers, chosen per calculation.

| Tier | Value error (relative to K) | Derivative error (relative to \|w′\|) |
|---|---|---|
| Fast (default) | ≤ 1 × 10⁻⁶ | ≤ 1 × 10⁻⁵ |
| Accurate | ≤ 2 × 10⁻¹³ | ≤ 2 × 10⁻¹² |

The bounds hold for y ∈ [0, 10⁶] and |x| ≤ 10⁸; points whose true K is below the smallest normal
double are exempt. Negative x follows from w(−x + iy) = conj(w(x + iy)) and
w′(−x + iy) = −conj(w′(x + iy)), applied exactly.

## Accurate tier

The value is a Rust port of Steven G. Johnson's Faddeeva Package (MIT licence; see
{doc}`/references`). The derivative is absrbrs's own, because the identity
w′ = −2zw + 2i/√π loses about |z|² in relative accuracy:

| Region | Derivative |
|---|---|
| \|z\| < 6 | the identity |
| x ≥ 6 and y ≤ 0.1 | derivative of exp(−z²) + continued fraction of depth 16 |
| \|z\| ≥ 6 otherwise | derivative of the continued fraction of depth 20 |

## Fast tier

A region-switched composite, chosen by measurement in phase 0
({doc}`/validation/faddeeva-prototype`). The first matching row wins:

| Region | Method |
|---|---|
| y ≤ 0.01 and 6 ≤ \|x\| < 27 | exp(−z²) + continued fraction, depth 5 |
| \|z\| ≥ 1750 | continued fraction, depth 0 (Lorentzian) |
| 50 ≤ \|z\| < 1750 | continued fraction, depth 1 (Humlíček region I) |
| 15 ≤ \|z\| < 50 | continued fraction, depth 2 |
| 9 ≤ \|z\| < 15 | continued fraction, depth 3 |
| 8 ≤ \|z\| < 9 | continued fraction, depth 4 |
| y ≤ 0.01 and \|x\| < 6 | Taylor series in iy, order 3, with Dawson's integral |
| otherwise | Weideman's rational approximation, 24 terms |

Each method returns the derivative of its own formula, so values and derivatives are mutually
consistent. In the Taylor region, which is not analytic in z, K_x and K_y are the formula's
partial derivatives.
```

Create `docs/validation/faddeeva-accuracy.md`:
````markdown
# Accuracy of the Faddeeva kernels

Every CI run checks both tiers against Arb values computed during the run (`validation/`): about
140 000 grid points, which include both sides of every Fast-tier boundary along and near the
real axis, and 200 000 random points over the validated domain. The errors use the metrics of
{doc}`/theory/faddeeva`. The Rust tests also check both tiers against the committed Arb grid,
and the Fast tier against the Accurate tier on 20 000 random points and around every region
boundary.

The latest measurement, written by `uv run --directory validation python -m absrbrs_validation
faddeeva`:

```{include} faddeeva-accuracy-table.md
```

The full report, with where each maximum occurs, is
{download}`faddeeva-accuracy.json <data/faddeeva-accuracy.json>`.
````

In `docs/index.md`:
- Replace the note's first sentence with "absrbrs is in early development: phase 1a (the Faddeeva kernels and `absrbrs.special`) is complete, and there is still no cross-section API."
- Add a toctree after "User guide":
  ````markdown
  ```{toctree}
  :maxdepth: 1
  :caption: Theory

  theory/faddeeva
  ```
  ````
- Add `validation/faddeeva-accuracy` as the first entry of the "Validation" toctree.

Change the README's status sentence to say phase 1a (Faddeeva kernels) is complete and phase 1b (line data, partition functions and the line-state layer) is next.

- [ ] **Step 6: Finish the follow-ups and the changelog**

Append to `docs/plans/phase-1a-followups.md`:
```markdown
## Accuracy (validation/, Arb computed at run time)

<copy docs/validation/faddeeva-accuracy-table.md here>

## For phase 1c (profile layer)

- `kernel_into` evaluates point by point. The profile layer evaluates one line's window, with
  sorted x and one y, so it can find each region's x-range once and run branch-free loops over
  them. Decide on that optimisation with the line benchmarks above.
- The Fast tier's margins are tightest at the depth-0 edge (|z| = 1750). If a future change to
  the depth-0 form costs accuracy, move the edge to 1800 (Lorentzian error ≈ 1.5/|z|²).
```
Replace `<copy … here>` with the actual table.

In `CHANGELOG.md` under `[Unreleased]` → `### Added`, add:
- "`absrbrs::faddeeva`: w(z) and its derivative, and the Voigt kernel triplet, at the Fast and Accurate tiers (the Accurate tier ports the MIT-licensed Faddeeva Package)."
- "`absrbrs.special.faddeeva` and `absrbrs.special.voigt_kernel`."
- "Validation of both tiers against Arb in CI; benchmarks against sasktran2's J-pole approximation."

- [ ] **Step 7: Run every check**

Run the whole "Checks" block of `CONTRIBUTING.md` and `uv run --locked pre-commit run --all-files`. Expected: everything passes, and the Sphinx build has no warnings.

- [ ] **Step 8: Commit**

```bash
git add validation .gitignore pyproject.toml .github CONTRIBUTING.md docs README.md CHANGELOG.md
git commit -m "feat(validation): both Faddeeva tiers against Arb at CI time; theory and accuracy docs

Co-Authored-By: <your model> <noreply@anthropic.com>"
```

---

## Phase 1a exit checklist (spec §13)

- [ ] Faddeeva kernels at both tiers, with derivatives: Tasks 1–3.
- [ ] Accuracy checked in CI against the committed Arb grid (Rust) and against dense Arb maps computed at run time (`validation`): Tasks 2, 3, 6.
- [ ] Faddeeva Package licence conditions met: notice in every ported file, `absrbrs/LICENSES/Faddeeva-MIT.txt` in `license-files`. Task 1.
- [ ] Python access through `absrbrs.special`: Task 5.
- [ ] Benchmarks against the J-pole recorded: Task 4.
- [ ] Docs: theory, accuracy, API. `docs/plans/phase-1a-followups.md` holds the decisions and the items for phase 1c.
