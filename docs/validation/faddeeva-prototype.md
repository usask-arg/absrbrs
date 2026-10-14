# Choosing the Fast-tier Faddeeva algorithm

Every Voigt line shape in absrbrs is the real part of the Faddeeva function
w(z) = exp(−z²) erfc(−iz), at z = x + iy with x = (ν − ν₀′)/σ_D and y = γ_L/σ_D. Here ν₀′ is the
pressure-shifted line centre, σ_D the Doppler 1/e half width and γ_L the Lorentz half width at
half maximum.

absrbrs will offer a Fast tier, with the targets below, and an Accurate tier (about 2 × 10⁻¹³, a
port of the Faddeeva Package). The Fast tier must hold the error of K = Re w to 1 × 10⁻⁶
relative, and that of each part of w′(z) to 1 × 10⁻⁵ relative to |w′|, over the whole domain: y
from 0 to 10⁶ and |x| up to 10⁸. Points where the true K is below the smallest normal double are
exempt.

Phase 0 chose the algorithm by measurement. The prototype in `prototypes/faddeeva/` evaluates
candidate approximations, each returning w and the derivative of its own formula:

- Weideman's rational approximation with 24, 32 or 40 terms;
- the Laplace continued fraction truncated at depth 0 to 16 (depth 0 is the Lorentzian and
  depth 1 Humlíček's region-I formula);
- exp(−z²) plus the continued fraction, for small y and moderate |x|;
- a Taylor series in y about the real axis, using Dawson's integral.

Four baselines are measured alongside:

- sasktran2's J-pole approximation (8 poles);
- HAPI's hum1wei24, re-implemented from the prototype's own methods: Humlíček's region I where
  |x| + y ≥ 15 and Weideman with 24 terms elsewhere;
- Weideman with 32 terms, and with 40 terms, each used alone.

The prototype measures them against Arb at about 164 000 points: the committed reference grid
(`reference/data/faddeeva.jsonl`) and a denser map computed with Arb, which samples both sides
of every region boundary the search can choose, along and near the real axis, where the errors
peak. A search over region-switched composites then picks the cheapest one that meets half of
each target, leaving room for the Rust implementation's rounding, and that also meets it on a
holdout of 400 000 random Arb points. Cost is a flop-count proxy, weighted by where the grid
points of representative atmospheric lines fall in the (|x|, y) plane: O₂ A-band and CO₂ 1.6 µm
lines at the surface, near 30 km and near 65 km, on a 0.005 cm⁻¹ grid with a 25 cm⁻¹ window.
Most of those points lie far from the line centre, where a shallow continued fraction is
accurate and cheap.

## Result

The selected composite uses the first region that contains the point; |x| is used throughout,
with w(−x + iy) = conj(w(x + iy)) and w′(−x + iy) = −conj(w′(x + iy)).

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

In the tables below, `cf n` is the continued fraction of depth n; `gauss_cf n` is exp(−z²) plus a
depth-n continued fraction; `taylor n` is the Taylor series in iy of order n, using Dawson's
integral; and `weideman n` is Weideman's rational approximation with n terms. "Expected cost" is
in flop-equivalents per point (a proxy), and "Workload share" is the fraction of the
representative lines' grid points that fall in the region.

These are the errors of the Python prototype; phase 1a validates the Rust kernel separately. The
derivative maximum of 1.8e-6 comes from the truncated Taylor series: its derivative is that of
the formula, and its y-part is a series one order shorter than the value's.

The maximum errors below are relative, as defined above. They are enormous for the baselines
because, far into the wings at zero or very low pressure, K is tiny (exp(−x²) or y/(√π x²)) and
a fixed-accuracy approximation's absolute error swamps it. HAPI's hum1wei24 also misses the
derivative target.

```{include} faddeeva-prototype-table.md
```

The full report, with every region's bounds, is
{download}`faddeeva_prototype.json <data/faddeeva_prototype.json>`. It is strict JSON, so `null`
stands for an unbounded region limit, or for an error that is infinite (a NaN somewhere).
Regenerate both files with

```bash
uv run --directory prototypes/faddeeva python -m faddeeva_prototype report
```

The dense-map and holdout inputs are not bit-portable across platforms, and the first run
computes about 555 000 Arb points (they are cached afterwards).

References are listed on the {doc}`/references` page.
