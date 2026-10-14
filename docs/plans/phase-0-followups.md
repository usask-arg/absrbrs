# Phase 0 follow-ups and decisions

Decisions taken in phase 0, and items deferred to later phases. This file is the most
up-to-date record. The spec (`docs/specs/2026-10-14-absrbrs-design.md`) is not edited to
match, apart from the error-type amendments to §7.4 and §9.2; the note under its title block
points here.

Phase 0 is complete at the merge of branch `phase-0-foundations`.

## Fast-tier Faddeeva algorithm (spec §6.3) — decided

The prototype (`prototypes/faddeeva/`, report in
`docs/validation/data/faddeeva_prototype.json`) selected this composite. Each point uses the
first region that contains it; |x| is used throughout, with w(−x + iy) = conj(w(x + iy)) and
w′(−x + iy) = −conj(w′(x + iy)).

| Region | Method |
|---|---|
| y ≤ 0.01 and 6 ≤ \|x\| < 27 | exp(−z²) + continued fraction, depth 5 |
| 8 ≤ \|z\| < 9 | continued fraction, depth 4 |
| 9 ≤ \|z\| < 15 | continued fraction, depth 3 |
| 15 ≤ \|z\| < 50 | continued fraction, depth 2 |
| 50 ≤ \|z\| < 1750 | continued fraction, depth 1 (Humlíček region I) |
| \|z\| ≥ 1750 | continued fraction, depth 0 (Lorentzian) |
| y ≤ 0.01 and \|x\| < 6 | Taylor series in iy, order 3 (needs Dawson's integral) |
| otherwise | Weideman, 24 terms |

The normative definitions of the methods are the prototype's functions in
`prototypes/faddeeva/faddeeva_prototype/candidates.py`: `gauss_cf`, `laplace_cf` (depth n:
coefficients k/2 for k = n, …, 1), `small_y_taylor` (including its derivative, which is that
of the truncated series in both x and y) and `weideman` (L = √(N/√2), coefficients by FFT).

Measured over 164 384 Arb grid points, which sample every region boundary along and near the
real axis where the errors peak, and confirmed on 400 000 random Arb points: maximum value
error 4.9e-7 and maximum derivative error 1.8e-6, against targets of 1e-6 and 1e-5. (The Task 5
review also checked full arcs at every |z| edge, y down to 1e-320, and fine x steps near the
axis.) The derivative maximum is in the Taylor region, where the derivative is that of the
formula: its y-part is a series one order shorter than the value's.

The expected cost proxy is 55 flop-equivalents per point, against 240 for sasktran2's J-pole
and 66 for HAPI's hum1wei24. About 98% of the points of atmospheric line windows fall in the
depth-0 and depth-1 continued-fraction regions, so the depth-0/depth-1 edge (1750) is the main
cost lever; the other choices are near-ties on the proxy (the cheapest valid composites differ
by under 1%).

Baselines, for the record: the J-pole (sasktran2) and hum1wei24 (HAPI) do not meet the value
target. Their relative errors explode where K is tiny, far in the wings at low pressure. This
confirms spec Appendix A, item 2.

## For phase 1a (Faddeeva kernels)

- Implement the composite above in Rust, returning K, K_x and K_y. In the continued-fraction
  regions (98% of points), evaluate the real and imaginary parts in real arithmetic, without
  complex division. Measure the real cost against the proxy: the proxy prices a complex
  division at 12 and a continued-fraction level at 40 flop-equivalents, so it ranks the
  composites but does not predict the speed-up.
- For a line with fixed y the regions become ranges of |x|: |z| ≥ r is |x| ≥ √(r² − y²), and
  the y ≤ 0.01 regions disappear when y > 0.01 (spec §6.4). Test the Rust accuracy map with
  the Rust classification itself: sqrt and hypot can round differently at an edge.
- The value margins are tightest at the depth-0 edge (4.9e-7 at |z| = 1750; the Lorentzian's
  relative error is about 1.5/|z|²), in the 15 ≤ |z| < 17.5 band (4.8e-7) and in the
  50 ≤ |z| < 60 band (4.0e-7), all against the 5e-7 half-target. Sample them densely at their
  edges. If rounding takes the depth-0 edge over 5e-7, move it to 1800; elsewhere add one level
  of depth.
- Taylor order 4 is a near-tie (about +0.02 on the proxy) and lowers the maximum derivative
  error from 1.8e-6 to about 5e-7; prefer it if the Rust cost confirms the tie.
- Dawson's integral D(x), for |x| < 6, is needed by every valid composite (no template point
  without the Taylor region passes). Its accuracy matters: the Taylor value term
  y (2/√π)(2x D(x) − 1) cancels as |x| → 6 and amplifies an error in D about 70 times. Measured
  in the prototype, D to 1e-8 relative gives 7e-7 (over the half-target), and D to 1e-9 gives
  8e-8. Use D to about 1e-9 relative or better: from the Accurate-tier port of the Faddeeva
  Package (its `w_im`), or compute 2x D(x) − 1 directly, without the cancellation.
- In the Taylor region, test K_y as ∂K/∂y of the formula, not as −Im of the x-derivative. The
  Taylor series uses w′(x) = −2xw + 2i/√π on the real axis, which amplifies relative error by
  about 2x² ≤ 72 for |x| < 6; that is acceptable with D to 1e-9, and it is why spec §6.1's rule
  against the identity does not apply there.
- `scipy.special.wofz` is Johnson's Faddeeva Package. It is a convenient independent check
  of the Accurate-tier Rust port, in addition to the Arb data.
- The dense map's and the holdout's inputs come from the platform's `10.0 **` and numpy's
  generator, and are not bit-portable. That does not matter for uncommitted maps, but any
  committed accuracy map must use portable inputs, as `reference/` does.

## Review notes for later phases

Recorded from the phase 0 task reviews.

- **I/O errors (phase 1b).** `Error::Io` keeps `kind` and a message, so Python gets
  `FileNotFoundError`/`PermissionError`/`OSError`, but `.errno` and `.filename` are `None`.
  When the file readers land, decide whether to also keep `raw_os_error: Option<i32>` and build
  `OSError(errno, strerror, path)` (standard text, CPython picks the subclass), or document that
  the attributes are unset.
- **uv cache-keys (phase 1b).** Data embedded in the core crate (`include_bytes!` of the TIPS
  tables) must be added to `[tool.uv] cache-keys`, or `uv run pytest` tests a stale extension;
  `rust-toolchain.toml` is missing from the keys as well.
- **Build profiles.** `lto = "fat"` and `codegen-units = 1` apply only to builds of this
  workspace; when sasktran2 depends on absrbrs, sasktran2's profiles apply. State this next to
  the spec §10 performance targets, and benchmark under both.
- **Free-threaded Python.** abi3 wheels exclude free-threaded builds (3.13t/3.14t); say so in
  the spec §9.3 packaging section.
- **Crate docs.** The crate-level doc in `absrbrs/src/lib.rs` points at a repository-relative
  spec path, which docs.rs cannot link; use the GitHub URL once the remote exists, and check the
  README links before the first crates.io publish.
- **Error variants.** Struct-variant fields are a semver commitment (enum-level
  `#[non_exhaustive]` does not cover them); acceptable at 0.0.x, revisit before 0.1.0.
- **Wheel notices for linked Rust dependencies (before the first wheel release, phase 1d).**
  Wheels statically link pyo3, pyo3-ffi, libc, once_cell and thiserror (MIT OR Apache-2.0).
  Add a cargo-deny licence allowlist and a generated notices file listed in `license-files`,
  and a row or paragraph in `THIRD_PARTY.md`.
- **Scheduled jobs (phase 1).** Spec §12 lists ignored (slow) tests and cross-code comparisons
  for `scheduled.yml`; add them when they exist.
- **Public paths.** The plans and spec mention local checkouts (`~/dev/scat-rs`,
  `~/dev/sasktran2`); refer to the GitHub repositories instead before publishing.

## Licence gates

Outcomes of the phase 0 licence checks (Task 6). The details and sources are in
`THIRD_PARTY.md`. Every row there is pending maintainer confirmation. Phase 1a may not merge
the Faddeeva Package port, and phase 1b may not commit the TIPS-2021 tables, until a maintainer
has confirmed those rows (name and date) in `THIRD_PARTY.md`.

- **Faddeeva Package (phase 1a): conditions.** It is MIT-licensed, so the Accurate-tier port
  may proceed once a maintainer has confirmed the row. The ported Rust file must carry the "Copyright (c) 2012 Massachusetts Institute
  of Technology" notice and the MIT permission notice in its header. Wheels contain no Rust
  source, and `pyproject.toml` has `license-files = ["LICENSE"]`. So in the same commit as the
  port, add a notice file inside the crate package (e.g. `absrbrs/LICENSES/Faddeeva-MIT.txt`)
  and list it in `pyproject.toml` `license-files`. Port from `Faddeeva.cc` only (SHA-256
  recorded in `THIRD_PARTY.md`), never from the ACM code for TOMS Algorithm 916 or 680
  (Poppe–Wijers).
- **TIPS-2021 (phase 1b): conditions.** The Zenodo record (10.5281/zenodo.4708099) is
  MIT-licensed, the paper (section 4) says the codes all use the MIT licence, and HAPI already
  redistributes the tables under MIT. The built-in tables may therefore be committed once a
  maintainer has confirmed the row. The
  HAPI-based builder is not required; user-supplied tables remain in scope (spec §5.5).
  Obligations:
  - Ship the copyright notice (Robert Gamache, 2021) and the MIT permission notice with the
    embedded data. In the same commit as the data, add a notice file inside the crate package
    (e.g. `absrbrs/LICENSES/TIPS-2021-MIT.txt`) and list it in `pyproject.toml`
    `license-files`, so the notice reaches both the crate and the wheel.
  - Cite Gamache et al. (2021).
  - Generate the tables with a script, recording for each file used Zenodo's published MD5
    and a locally computed SHA-256.

  Primary source: the `.QTpy` files. They carry 6–8 significant figures (H2O isotopologue 1:
  `'1.32837'` at 10 K, `'174.58130'` at 296 K), and HAPI's embedded tables carry about 7. The
  Fortran DATA tables in `TIPS_2021.zip` (`BD_TIPS_2021.FOR`) carry only 5 (`1.3284E+00`).
  Building from the Fortran would add up to about 5e-5 relative rounding error, and could not
  meet spec §5.5's "equal HAPI's partitionSum exactly at nodes".

  The `.QTpy` files are third-party Python pickles, and a matching hash does not make
  unpickling them safe. Inspect each one with `pickletools` first. `1_1.QTpy` is a plain dict
  of str to str values, with no `GLOBAL` or `REDUCE` opcodes. Then load it, offline, with a
  restricted `pickle.Unpickler` whose `find_class` refuses all globals. Never unpickle anything
  at runtime. Use the Fortran tables as an independent cross-check: the two must agree to the
  Fortran's 5 significant figures.
- **TIPS grids and HAPI version (phase 1b).** HAPI 1.3.0.0's `partitionSum` defaults to
  `version=2025` (TIPS-2025). Its embedded TIPS-2021 grid is 1, 10, 20, … K, whereas the Zenodo
  tables use 1 K steps up to 20 K and 2 K steps above. Spec §5.5's "tabulated at 1 K steps" and
  "equal HAPI's partitionSum exactly at nodes" need revisiting in phase 1b. Pin `version=2021`
  and compare at HAPI's nodes.
- **HITRAN line data (phase 1b): blocked.** HITRAN states no licence or terms of use. It is
  free of charge and asks for citation, but does not address redistribution, so it counts as
  all rights reserved.
  - Do not commit, ship or publish HITRAN records, including in CI artifacts and caches. This
    does not cover users, or the spec §8.5 benchmark jobs, downloading HITRAN data at runtime.
  - The loader golden tests use synthetic `.par` records instead of HITRAN snippets (this
    replaces spec §8.6). The synthetic records must exercise the same fields and edge cases.
  - To unblock: an explicit licence or written permission from the HITRAN team.
- **HAPI: cleared** as an optional, user-installed dependency of the table builder (if built)
  and as a development comparison. No HAPI code is copied.
- **sasktran2 J-pole constants: conditions.** They are used only by the prototype, which is in
  the repository but not in the crate or wheel. The comment above `_JPOLE_B` in
  `candidates.py` carries the provenance and copyright notice; keep it.
- **Weideman (1994): cleared.** The prototype follows the published algorithm, whose few-line
  listing is the direct expression of the formula; low risk. Phase 1a derives the Rust from the
  paper's equations. Cite the paper.
- **python-flint and FLINT: cleared.** They are development tools only; FLINT is
  LGPL-3.0-or-later and must never be linked or shipped.
- **Aside.** TIPS-2025 (for HITRAN2024; Zenodo 10.5281/zenodo.17191976) is CC-BY-4.0 according
  to its Zenodo record. If phase 1b would rather use it than TIPS-2021, check it and register it
  first.
