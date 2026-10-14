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
            for (depth, (w, dw)) in [(0, lorentzian(x, y)), (1, humlicek_1(x, y))] {
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
