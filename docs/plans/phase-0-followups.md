# Phase 0 follow-ups and decisions

Decisions taken in phase 0, and items deferred to later phases. This file is the most
up-to-date record; the spec (`docs/specs/2026-10-14-absrbrs-design.md`) is not edited to
match.

## Fast-tier Faddeeva algorithm (spec §6.3) — decided

The prototype (`prototypes/faddeeva/`, report in
`docs/validation/data/faddeeva_prototype.json`) selected this composite. Each point uses the
first region that contains it; |x| is used throughout, with w(−x + iy) = conj(w(x + iy)).

| Region | Method |
|---|---|
| y ≤ 0.01 and 6 ≤ \|x\| < 27 | exp(−z²) + continued fraction, depth 5 |
| 8 ≤ \|z\| < 10 | continued fraction, depth 4 |
| 10 ≤ \|z\| < 15 | continued fraction, depth 3 |
| 15 ≤ \|z\| < 50 | continued fraction, depth 2 |
| 50 ≤ \|z\| < 3000 | continued fraction, depth 1 (Humlíček region I) |
| \|z\| ≥ 3000 | continued fraction, depth 0 (Lorentzian) |
| y ≤ 0.01 and \|x\| < 6 | Taylor series in iy, order 3 (needs Dawson's integral) |
| otherwise | Weideman, 24 terms |

Measured over 156 802 Arb grid points, including every region boundary, and confirmed on
400 000 random Arb points: maximum value error 4.8e-7 and maximum derivative error 1.8e-6,
against targets of 1e-6 and 1e-5. The derivative maximum is in the Taylor region, where the
derivative is that of the formula: its y-part is a series one order shorter than the value's. The expected cost proxy is 61 flop-equivalents
per point, against 240 for sasktran2's J-pole and 66 for HAPI's hum1wei24. About 98% of the
points of atmospheric line windows fall in the depth-0 and depth-1 continued-fraction regions.

Baselines, for the record: the J-pole (sasktran2) and hum1wei24 (HAPI) do not meet the value
target. Their relative errors explode where K is tiny, far in the wings at low pressure. This
confirms spec Appendix A, item 2.

## For phase 1a (Faddeeva kernels)

- Implement the composite above in Rust, returning K, K_x and K_y. In the continued-fraction
  regions (98% of points), evaluate the real and imaginary parts in real arithmetic, without
  complex division. Measure the real cost against the proxy.
- The value margins are tightest in the 15 ≤ |z| < 20 band (4.8e-7 against the 5e-7
  half-target) and the 50 ≤ |z| < 100 band (4.0e-7, at |z| = 50 and small y). The Rust accuracy
  map must sample both densely, at their edges especially. If rounding takes either over 1e-6,
  add one level of depth there.
- The Taylor region needs Dawson's integral D(x) for |x| < 6. Take it from the
  Accurate-tier port of the Faddeeva Package (its `w_im`), or use a dedicated rational
  approximation. Decide when the Accurate tier is ported.
- `scipy.special.wofz` is Johnson's Faddeeva Package. It is a convenient independent check
  of the Accurate-tier Rust port, in addition to the Arb data.
- The dense map's inputs come from the platform's `10.0 **` and are not portable. That does
  not matter for an uncommitted map, but any committed accuracy map must use portable inputs,
  as `reference/` does.
