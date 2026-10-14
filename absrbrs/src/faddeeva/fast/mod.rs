//! The Fast tier (spec §6.3): the composite chosen in phase 0 (docs/plans/phase-0-followups.md).

mod taylor;
mod weideman;

use num_complex::Complex64;

use super::FRAC_1_SQRT_PI;
use super::cf::{gauss_cf, laplace_cf};

/// Near-axis threshold: y ≤ Y_NEAR_AXIS selects the Gauss continued fraction or the Taylor series.
pub(crate) const Y_NEAR_AXIS: f64 = 0.01;
/// Near the axis, the Taylor series for x < X_TAYLOR_GAUSS, the Gauss continued fraction above.
pub(crate) const X_TAYLOR_GAUSS: f64 = 6.0;
/// Near the axis, the Gauss continued fraction for x < X_GAUSS_MAX.
pub(crate) const X_GAUSS_MAX: f64 = 27.0;
/// Laplace continued fraction of depth 4 for R_CF4 ≤ |z| < R_CF3.
pub(crate) const R_CF4: f64 = 8.0;
/// Depth 3 for R_CF3 ≤ |z| < R_CF2.
pub(crate) const R_CF3: f64 = 9.0;
/// Depth 2 for R_CF2 ≤ |z| < R_HUMLICEK1.
pub(crate) const R_CF2: f64 = 15.0;
/// Depth 1 (Humlíček's region I) for R_HUMLICEK1 ≤ |z| < R_LORENTZIAN.
pub(crate) const R_HUMLICEK1: f64 = 50.0;
/// Depth 0 (Lorentzian) for |z| ≥ R_LORENTZIAN.
pub(crate) const R_LORENTZIAN: f64 = 1750.0;

/// The regions of the composite, each named after its method.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Region {
    /// y ≤ Y_NEAR_AXIS and X_TAYLOR_GAUSS ≤ x < X_GAUSS_MAX: `gauss_cf(z, 5)`.
    GaussCf5,
    /// |z| ≥ R_LORENTZIAN: depth 0, in real arithmetic.
    Lorentzian,
    /// R_HUMLICEK1 ≤ |z| < R_LORENTZIAN: depth 1, in real arithmetic.
    Humlicek1,
    /// R_CF2 ≤ |z| < R_HUMLICEK1: `laplace_cf(z, 2)`.
    Cf2,
    /// R_CF3 ≤ |z| < R_CF2: `laplace_cf(z, 3)`.
    Cf3,
    /// R_CF4 ≤ |z| < R_CF3: `laplace_cf(z, 4)`.
    Cf4,
    /// y ≤ Y_NEAR_AXIS and x < X_TAYLOR_GAUSS (and |z| < R_CF4): the order-3 Taylor series.
    Taylor3,
    /// Everything else: Weideman's 24-term rational approximation.
    Weideman24,
}

/// The region of the composite at x + iy, for x ≥ 0 and y ≥ 0 (the first matching row of the
/// table in docs/plans/phase-0-followups.md wins), and r = |z| = `x.hypot(y)`, which the
/// Lorentzian reuses (the compiler does not merge two `hypot` calls).
#[inline]
pub(crate) fn region(x: f64, y: f64) -> (Region, f64) {
    let r = x.hypot(y);
    let region = if y <= Y_NEAR_AXIS && (X_TAYLOR_GAUSS..X_GAUSS_MAX).contains(&x) {
        Region::GaussCf5
    } else if r >= R_LORENTZIAN {
        Region::Lorentzian
    } else if r >= R_HUMLICEK1 {
        Region::Humlicek1
    } else if r >= R_CF2 {
        Region::Cf2
    } else if r >= R_CF3 {
        Region::Cf3
    } else if r >= R_CF4 {
        Region::Cf4
    } else if y <= Y_NEAR_AXIS {
        Region::Taylor3
    } else {
        Region::Weideman24
    };
    (region, r)
}

/// w(z) and w′(z) for Re z ≥ 0 and Im z ≥ 0, by the method of [`region`]. In the Taylor region
/// the series is not analytic in z, and the "w′" returned is K_x − iK_y of that formula: its
/// imaginary part is −K_y, not ∂ Im w/∂x.
pub(crate) fn w_and_derivative(z: Complex64) -> (Complex64, Complex64) {
    let (x, y) = (z.re, z.im);
    let (region, r) = region(x, y);
    match region {
        Region::GaussCf5 => gauss_cf(z, 5),
        Region::Lorentzian => lorentzian(x, y, r),
        Region::Humlicek1 => humlicek_1(x, y),
        Region::Cf2 => laplace_cf(z, 2),
        Region::Cf3 => laplace_cf(z, 3),
        Region::Cf4 => laplace_cf(z, 4),
        Region::Taylor3 => taylor::w_and_derivative(x, y),
        Region::Weideman24 => weideman::w_and_derivative(z),
    }
}

/// Depth 0: w = i/(√π z), w′ = −i/(√π z²). With r = |z|, c = x/r and s = y/r:
/// w = (s + ic)/(√π r) and w′ = −(2cs + i(c − s)(c + s))/(√π r²), dividing by r twice so that
/// nothing overflows for |z| up to f64::MAX.
fn lorentzian(x: f64, y: f64, r: f64) -> (Complex64, Complex64) {
    let (c, s) = (x / r, y / r);
    let k = FRAC_1_SQRT_PI / r;
    let k2 = k / r;
    (
        Complex64::new(s * k, c * k),
        Complex64::new(-2.0 * c * s * k2, -(c - s) * (c + s) * k2),
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

    /// The next double above a positive v.
    fn up(v: f64) -> f64 {
        f64::from_bits(v.to_bits() + 1)
    }

    /// The next double below a positive v.
    fn down(v: f64) -> f64 {
        f64::from_bits(v.to_bits() - 1)
    }

    #[test]
    fn region_follows_the_table() {
        use Region::*;
        let interior = [
            ((10.0, 0.001), GaussCf5),
            ((5000.0, 3.0), Lorentzian),
            ((100.0, 20.0), Humlicek1),
            ((20.0, 5.0), Cf2),
            ((10.0, 5.0), Cf3),
            ((6.0, 5.5), Cf4),
            ((2.0, 0.001), Taylor3),
            ((2.0, 1.0), Weideman24),
        ];
        // (below the edge, at the edge): points on both sides of every boundary.
        let edges = [
            // The |z| circles, off the near-axis strip.
            (
                (0.0, down(R_LORENTZIAN)),
                Humlicek1,
                (0.0, R_LORENTZIAN),
                Lorentzian,
            ),
            ((0.0, down(R_HUMLICEK1)), Cf2, (0.0, R_HUMLICEK1), Humlicek1),
            ((0.0, down(R_CF2)), Cf3, (0.0, R_CF2), Cf2),
            ((0.0, down(R_CF3)), Cf4, (0.0, R_CF3), Cf3),
            ((0.0, down(R_CF4)), Weideman24, (0.0, R_CF4), Cf4),
            // On the near-axis strip, beyond X_GAUSS_MAX.
            ((down(R_HUMLICEK1), 0.0), Cf2, (R_HUMLICEK1, 0.0), Humlicek1),
            (
                (down(R_LORENTZIAN), 0.0),
                Humlicek1,
                (R_LORENTZIAN, 0.0),
                Lorentzian,
            ),
            // The x thresholds of the strip.
            (
                (down(X_TAYLOR_GAUSS), 0.0),
                Taylor3,
                (X_TAYLOR_GAUSS, 0.0),
                GaussCf5,
            ),
            ((down(X_GAUSS_MAX), 0.0), GaussCf5, (X_GAUSS_MAX, 0.0), Cf2),
            // The y threshold: (at the edge, above it) here, so the order is reversed.
            (
                (2.0, Y_NEAR_AXIS),
                Taylor3,
                (2.0, up(Y_NEAR_AXIS)),
                Weideman24,
            ),
            ((10.0, Y_NEAR_AXIS), GaussCf5, (10.0, up(Y_NEAR_AXIS)), Cf3),
            (
                (7.0, Y_NEAR_AXIS),
                GaussCf5,
                (7.0, up(Y_NEAR_AXIS)),
                Weideman24,
            ),
        ];
        for ((x, y), want) in interior {
            assert_eq!(region(x, y).0, want, "({x}, {y})");
        }
        for ((x0, y0), want0, (x1, y1), want1) in edges {
            assert_eq!(region(x0, y0).0, want0, "({x0}, {y0})");
            assert_eq!(region(x1, y1).0, want1, "({x1}, {y1})");
        }
    }

    /// The Lorentzian does not overflow: w and w′ stay finite up to |z| = f64::MAX.
    #[test]
    fn lorentzian_is_finite_for_huge_z() {
        for (x, y) in [(1e200, 0.0), (1e155, 1e155), (0.0, 1e160), (f64::MAX, 0.0)] {
            let (w, dw) = w_and_derivative(Complex64::new(x, y));
            assert!(
                w.re.is_finite() && w.im.is_finite() && dw.re.is_finite() && dw.im.is_finite(),
                "({x}, {y}): w = {w}, w′ = {dw}"
            );
        }
    }

    /// The real-arithmetic forms equal the complex continued fraction of the same depth.
    #[test]
    fn real_forms_match_the_continued_fraction() {
        for (x, y) in [
            (50.0, 0.0),
            (60.0, 1e-3),
            (200.0, 37.0),
            (1749.0, 5.0),
            (1750.0, 0.0),
            (1e5, 1e4),
            (3e3, 2e3),
        ] {
            let z = Complex64::new(x, y);
            for (depth, (w, dw)) in [(0, lorentzian(x, y, x.hypot(y))), (1, humlicek_1(x, y))] {
                let (wc, dwc) = laplace_cf(z, depth);
                assert!(
                    (w.re - wc.re).abs() <= 1e-14 * wc.re.abs().max(f64::MIN_POSITIVE),
                    "K, depth {depth} at ({x}, {y})"
                );
                assert!(
                    (w - wc).norm() <= 1e-14 * wc.norm(),
                    "w, depth {depth} at ({x}, {y})"
                );
                assert!(
                    (dw - dwc).norm() <= 1e-14 * dwc.norm(),
                    "w′, depth {depth} at ({x}, {y})"
                );
            }
        }
    }
}
