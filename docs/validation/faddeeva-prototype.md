# Choosing the Fast-tier Faddeeva algorithm

Every Voigt line shape in absrbrs is the real part of the Faddeeva function
w(z) = exp(−z²) erfc(−iz), at z = x + iy with x = (ν − ν₀′)/σ_D and y = γ_L/σ_D. The Fast
accuracy tier must hold the error of K = Re w to 1 × 10⁻⁶ relative, and that of each part of
w′(z) to 1 × 10⁻⁵ relative to |w′|, over the whole domain: y from 0 to 10⁶ and |x| up to 10⁸.
Values below the smallest normal double are exempt.

Phase 0 chose the algorithm by measurement. The prototype in `prototypes/faddeeva/` evaluates
candidate approximations, each returning w and the derivative of its own formula:

- Weideman's rational approximation with 24, 32 or 40 terms;
- the Laplace continued fraction truncated at depth 0 to 16 (depth 0 is the Lorentzian and
  depth 1 Humlíček's region-I formula);
- exp(−z²) plus the continued fraction, for small y and moderate |x|;
- a Taylor series in y about the real axis, using Dawson's integral;
- sasktran2's 8-pole approximation, and HAPI's hum1wei24, as baselines.

It measures them against Arb at about 164 000 points: the committed reference grid
(`reference/data/faddeeva.jsonl`) and a denser map computed with Arb, which samples every region
boundary the search can choose. A search over region-switched composites then picks the
cheapest one that meets half of each target, leaving room for the Rust implementation's
rounding, and that also meets it on a holdout of 400 000 random Arb points. Cost is a flop-count proxy, weighted by where the
grid points of representative atmospheric lines fall in the (|x|, y) plane: O₂ A-band and CO₂
1.6 µm lines at the surface, near 30 km and near 65 km, on a 0.005 cm⁻¹ grid with a 25 cm⁻¹
window. Most of those points lie far from the line centre, where a shallow continued fraction
is accurate and cheap.

The maximum errors below are relative, as defined above. They are enormous for the baselines
because, far into the wings at zero or very low pressure, K is tiny (exp(−x²) or y/x²) and a
fixed-accuracy approximation's absolute error swamps it.

```{include} faddeeva-prototype-table.md
```

The full report, with every region's bounds, is
`docs/validation/data/faddeeva_prototype.json`. It is strict JSON, so `null` stands for an
unbounded region limit, or for an error that is infinite (a NaN somewhere). Regenerate both files with

```bash
uv run --directory prototypes/faddeeva python -m faddeeva_prototype report
```
