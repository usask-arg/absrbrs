# Phase 1a follow-ups and decisions

## Benchmarks (`cargo bench -p absrbrs --bench faddeeva`)

Measured on an Apple M4 Pro (arm64), macOS 26.7.1, 2026-10-14, release profile (fat LTO, one
codegen unit). Times are criterion's median estimates. The region benchmarks time one fixed
(x, y) per region, so the branch predictor sees the same path every iteration; the line-window
benchmarks are the realistic workload.

| Region | Fast (ns/point) | Accurate (ns/point) | J-pole (ns/point) |
|---|---|---|---|
| weideman | 48.46 | 35.77 | 14.61 |
| taylor | 15.25 | 34.65 | 14.53 |
| exp_cf5 | 22.37 | 98.95 | 14.75 |
| cf4 | 14.63 | 118.53 | 14.81 |
| cf3 | 11.60 | 111.54 | 15.11 |
| cf2 | 8.80 | 100.91 | 14.84 |
| cf1 | 4.94 | 95.30 | 14.78 |
| cf0 | 5.17 | 93.50 | 19.66 |

Line windows are 10 001 grid points each (0.005 cm⁻¹ over 25 cm⁻¹); the times are the median
per window divided by 10 001.

| Line window | Fast (ns/point) | Accurate (ns/point) | J-pole (ns/point) |
|---|---|---|---|
| o2_surface | 7.65 | 144.83 | 21.83 |
| o2_30km | 6.97 | 115.59 | 15.26 |
| o2_65km | 4.82 | 77.16 | 15.22 |
| co2_surface | 4.47 | 94.24 | 15.26 |
| co2_30km | 6.72 | 148.94 | 23.42 |
| co2_65km | 4.52 | 77.45 | 15.20 |

On line windows the Fast tier takes 0.29 to 0.46 of the J-pole's time per point, that is 2.2 to
3.5 times faster. The phase 0 cost proxy predicted about 4 times fewer flop-equivalents, so the
measured gain is somewhat smaller than predicted, but in the same range. The depth-0 and depth-1
real-arithmetic forms dominate as expected: cf0 and cf1 cost about 5 ns per point, the cheapest
of any region and about 3 times cheaper than the J-pole's fixed 15 ns, and the windows that run
closest to those forms (co2_surface, o2_65km, co2_65km, at about 4.5 ns) are the fastest. The
windows that spend more points in the cf2 to cf4, exp_cf5, taylor and Weideman regions (o2_30km,
co2_30km, o2_surface, at 6.7 to 7.7 ns) are slower, and the Weideman region on its own is about
3 times slower than the J-pole, but few grid points fall in it. The Accurate tier costs 77 to
149 ns per point, 16 to 22 times the Fast tier on these windows. The per-point `hypot` in the
Fast tier's classifier is significant: in a standalone loop on this machine `hypot` costs about
1.4 ns per point against 0.1 ns for x² + y², which is roughly a quarter to a third of the
cf0/cf1 point times, so a phase 1c option is to compare x² + y² against squared region edges and
call `hypot` only where the Lorentzian needs it.
