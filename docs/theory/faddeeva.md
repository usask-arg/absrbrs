# The Faddeeva function and the Voigt kernel

Every Voigt line shape is K(x, y) = Re w(x + iy), the real part of the Faddeeva function
w(z) = exp(−z²) erfc(−iz), at x = (ν − ν₀′)/σ_D and y = γ_L/σ_D ≥ 0 in the line-shape model.
The Jacobians need its partial derivatives: K_x = Re w′ and K_y = −Im w′, since w is analytic.
absrbrs computes the triplet (K, K_x, K_y) at two accuracy tiers, chosen per calculation.

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

So for |z| ≥ 6 the value (Johnson's) and the derivative come from different formulas, a
deliberate exception to spec §6.1's "same formula": the value is accurate to about 1e-13 and
the derivative to a few 1e-13 for |z| ≥ 6, so they are consistent to that level.

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
