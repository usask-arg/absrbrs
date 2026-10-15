# Phase 1a follow-ups and decisions

## Benchmarks (`cargo bench -p absrbrs --bench faddeeva`)

Measured on an Apple M4 Pro (arm64), macOS 26.7.1, 2026-10-14, release profile (fat LTO, one
codegen unit). Times are criterion's median estimates, and each cell is the lower of two clean
runs. The region benchmarks time one fixed (x, y) per region, so the branch predictor sees the
same path every iteration; the line-window benchmarks are the realistic workload.

Measurement protocol: the J-pole does no data-dependent work, so it must cost the same in every
row. Run the whole bench twice (`-- --save-baseline a`, then `b`) and accept a run only if
every J-pole row is flat (14.5 to 15.5 ns per point). Background load shows up as bursts in
some rows (J-pole rows at 20 to 23 ns, or Fast and Accurate rows up to 50% high); rerun until
two runs are clean. Here the first two runs were contaminated by another process on the machine
and the third and fourth were clean, with J-pole at 14.54 to 14.82 ns on the region benchmarks
and 14.98 to 15.25 ns on the line windows. The line closures end with
`black_box((&mut k, &mut kx, &mut ky))`, so the stores to the output arrays cannot be
eliminated.

| Region | Fast (ns/point) | Accurate (ns/point) | J-pole (ns/point) |
|---|---|---|---|
| weideman | 49.01 | 35.65 | 14.70 |
| taylor | 15.17 | 34.80 | 14.68 |
| exp_cf5 | 22.54 | 98.57 | 14.56 |
| cf4 | 14.60 | 116.01 | 14.61 |
| cf3 | 11.59 | 111.30 | 14.73 |
| cf2 | 8.77 | 101.09 | 14.70 |
| cf1 | 4.93 | 93.81 | 14.54 |
| cf0 | 5.16 | 89.86 | 14.58 |

Line windows are 10 001 grid points each (0.005 cm⁻¹ over 25 cm⁻¹); the times are the median
per window divided by 10 001.

| Line window | Fast (ns/point) | Accurate (ns/point) | J-pole (ns/point) |
|---|---|---|---|
| o2_surface | 4.84 | 91.10 | 14.98 |
| o2_30km | 4.83 | 75.85 | 15.03 |
| o2_65km | 4.79 | 74.93 | 15.01 |
| co2_surface | 4.42 | 92.31 | 15.11 |
| co2_30km | 4.50 | 92.91 | 15.15 |
| co2_65km | 4.48 | 75.94 | 15.13 |

The J-pole baseline is the 8-pole sum in plain complex arithmetic. sasktran2's own derivative
path precomputes per-line pole-pair coefficients and uses f64x8 SIMD under its nightly-only
feature, which this baseline does not represent (an informal scratch measurement, not
reproducible from the repository, put a scalar pole-paired version at 12.9 ns per point against
13.4 ns for this one), so the comparison is with scalar J-pole code.

All six line windows cost about the same on the Fast tier, 4.4 to 4.9 ns per point, because
more than 95% of their points use the depth-0 and depth-1 real forms. By occupancy, o2_30km and
o2_65km are 96 to 97% Humlicek depth 1; co2_surface is 44% depth 1 and 55% Lorentzian; the
Weideman region has at most 0.47% of the points in any window. That matches the point benchmarks:
cf0 and cf1 cost about 5 ns, the cheapest regions. The Fast tier is 3.1 to 3.4 times faster than
the J-pole on these windows (the phase 0 cost proxy predicted about 4 times fewer
flop-equivalents), and the Accurate tier costs 75 to 92 ns per point, 16 to 21 times the Fast
tier. Among the other regions the Fast tier is faster than the J-pole in cf2 and cf3, level with it
in cf4 and taylor, 1.5 times slower in exp_cf5 and 3.3 times slower in the Weideman region (see
below), but few points of these line windows fall there.

Per-point `hypot` in the Fast tier's classifier: in a standalone loop on this machine `hypot`
costs about 1.4 ns per point against 0.1 ns for x² + y². That is about a quarter of the cf0/cf1
time (4.9 to 5.2 ns), but it is an upper-bound estimate from a standalone loop and should be
confirmed in phase 1c. A phase 1c option is to compare x² + y² against squared region edges.
The Lorentzian reuses r, so a squared-edge classifier still needs r (or 1/r²) in that region, and
saves the `hypot` only where r is not used.

## Accuracy (validation/, Arb computed at run time)

| Tier | Points | Max value error | Target | Max derivative error | Target |
|---|---|---|---|---|---|
| fast | 347825 | 4.9e-07 | 1e-06 | 1.8e-06 | 1e-05 |
| accurate | 347825 | 9.6e-14 | 2e-13 | 8.8e-13 | 2e-12 |

The grid includes Johnson's own branch thresholds (|x| = 10, 28, 4000, 1e7; y = 1e-10, 5, 7) and
the double on the other side of each.

## Decisions

- **Accurate-tier derivative (spec §6.1 exception).** For |z| ≥ 6 the derivative comes from an
  asymptotic formula (exp(−z²) + continued fraction of depth 16 near the axis, continued fraction
  of depth 20 elsewhere) while the value is Johnson's. Both are accurate to about 1e-13, so value
  and derivative are consistent to that level; measured on 565k Arb points: value ≤ 9.7e-14,
  derivative ≤ 9.3e-13 (tightest: the identity just inside |z| = 6).
- **Johnson's large-x second loop** (`np = n0 + dn++`, then `exp(-sqr(a*dn+dx))`) has an index
  mismatch in the source, ported as written: it is unreachable on the machine-precision path.

## For phase 1c

- Weideman region speed. At (2, 0.5) the Fast tier takes 49 ns and the Accurate tier 36 ns.
  Across the Weideman region (y > 0.01, |z| < 8): for |z| < 2, Fast takes about 50 to 67 ns
  against about 43 ns for Accurate; for 2 to 6, about 50 against 45 to 47; for 6 to 8, about 50
  against 129; over the whole region about 50 against 83. So routing the whole region to
  Accurate would lose. Options: an Estrin or even/odd split Horner to break the 23-step
  dependency chain; a smaller N, if the 1e-6 and 1e-5 error bounds still hold (re-validate);
  or routing |z| ≲ 5 to Accurate. Occupancy is 0 to 0.47% in these windows (at most about 5% of
  line time), but it matters for densely sampled broad lines.
- `kernel_into` evaluates point by point. The profile layer evaluates one line's window, with
  sorted x and one y, so it can find each region's x-range once and run branch-free loops over
  them. Decide on that optimisation with the line benchmarks above.
- The Fast tier's margins are tightest at the depth-0 edge (|z| = 1750). If a future change to
  the depth-0 form costs accuracy, move the edge to 1800 (Lorentzian error ≈ 1.5/|z|²).

- Accurate tier cost (spec §10.2 sets "Accurate tier ≤ 4 × Fast tier"): measured 16 to 21 times on
  the line windows (75 to 92 against 4.4 to 4.9 ns per point). Likely cause: for |z| ≥ 6 the
  derivative always runs a continued fraction of depth 16 or 20 (about 20 complex inversions),
  even where Johnson's value needs ν ≤ 2. Option: choose the derivative's depth by |z| for 1e-13,
  as the Fast tier does for 1e-6, re-validate, and either meet the target or record a spec change.
- Accuracy maps: spec §8.2 and §12 call for published accuracy maps, but the accuracy page has
  overall maxima only. Add per-region maxima (binning by the Fast region classifier) or a heatmap
  in phase 1d.

## For phase 1b

The open items of the "Review notes" in `docs/plans/phase-0-followups.md`: I/O errors carrying the
errno and file name, the uv cache keys (including `rust-toolchain.toml` and the embedded data),
and the TIPS notice file and the safe `.QTpy` conversion, under the Licence gates.

## Phase 0 items for phase 1a: outcome

- Taylor order 4 (a near-tie with order 3 in phase 0): order 3 kept.
- The `scipy.special.wofz` cross-check: superseded by the validation against Arb; not done.
- Real-arithmetic forms of the continued fractions: implemented for depths 0 and 1 only; depths 2
  to 5 use complex arithmetic.
- Dawson accuracy: taken from the Accurate port's `w_im`, at most 4.3e-15.

## Before the first release

- Notices for statically linked Rust crates. The wheel statically links the Rust crates that
  `absrbrs-py` and `absrbrs` depend on: pyo3, rust-numpy (BSD-2-Clause), ndarray, matrixmultiply,
  num-integer, rawpointer, rustc-hash and portable-atomic-util (MIT or Apache-2.0), num-complex,
  thiserror, and the rest of the lockfile's runtime tree. Phase 1b adds serde_json, now a core
  dependency, with serde_core and itoa (MIT or Apache-2.0), memchr (Unlicense or MIT) and zmij
  (MIT). Their licences require reproducing the
  copyright and permission notices in binary redistributions. Before the first wheel is published,
  add a generated notice file for the bundled crates (for example with `cargo-about`), list it in
  `pyproject.toml` `license-files`, and add a paragraph to `THIRD_PARTY.md`; or record an explicit
  decision not to.
