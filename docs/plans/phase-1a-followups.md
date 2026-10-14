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
feature, which this baseline does not represent (a scalar pole-paired version measured 12.9 ns
per point against 13.4 ns for this one), so the comparison is with scalar J-pole code.

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

## For phase 1c

- Weideman region speed. At (2, 0.5) the Fast tier takes 49 ns and the Accurate tier 36 ns.
  Across the Weideman region (y > 0.01, |z| < 8): for |z| < 2, Fast takes about 50 to 67 ns
  against about 43 ns for Accurate; for 2 to 6, about 50 against 45 to 47; for 6 to 8, about 50
  against 129; over the whole region about 50 against 83. So routing the whole region to
  Accurate would lose. Options: an Estrin or even/odd split Horner to break the 23-step
  dependency chain; a smaller N, if the 1e-6 and 1e-5 error bounds still hold (re-validate);
  or routing |z| ≲ 5 to Accurate. Occupancy is 0 to 0.47% in these windows (at most about 5% of
  line time), but it matters for densely sampled broad lines.
