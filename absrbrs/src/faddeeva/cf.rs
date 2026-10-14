//! The Laplace continued fraction for w(z), and exp(−z²) plus it, each with the derivative of
//! its own formula (forward differentiation of the recurrence).

use num_complex::Complex64;

use super::FRAC_1_SQRT_PI;

const I: Complex64 = Complex64::new(0.0, 1.0);

/// w(z) ≈ (i/√π) / (z − ½/(z − 1/(z − (3/2)/(z − …)))) truncated after `depth` levels, and the
/// derivative of that truncation. Depth 0 is the Lorentzian i/(√π z); depth 1 is Humlíček's
/// region-I formula. An asymptotic expansion: for large |z| only (its truncations have poles on
/// the real axis within |x| < 6 for depth ≤ 20).
pub(crate) fn laplace_cf(z: Complex64, depth: u32) -> (Complex64, Complex64) {
    let mut f = z;
    let mut fp = Complex64::new(1.0, 0.0);
    for k in (1..=depth).rev() {
        let a = 0.5 * f64::from(k);
        let inv = f.inv();
        (f, fp) = (z - a * inv, 1.0 + a * fp * inv * inv);
    }
    let inv = f.inv();
    (
        I * FRAC_1_SQRT_PI * inv,
        -I * FRAC_1_SQRT_PI * fp * inv * inv,
    )
}

/// exp(−z²) plus the continued fraction, and its derivative: keeps the Gaussian core that the
/// asymptotic expansion lacks. exp(−z²) overflows for y above about 26.6: small y only.
pub(crate) fn gauss_cf(z: Complex64, depth: u32) -> (Complex64, Complex64) {
    let g = (-(z * z)).exp();
    let (w, dw) = laplace_cf(z, depth);
    // z·g first: for |x| ≥ 2^1023, 2z overflows while g underflows to 0, and (2z)·g is ∞·0 = NaN.
    (w + g, dw - 2.0 * (z * g))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shallow_fractions_are_the_classical_formulas() {
        let z = Complex64::new(3.0, 2.0);
        let (w0, dw0) = laplace_cf(z, 0);
        assert!((w0 - I * FRAC_1_SQRT_PI / z).norm() <= 1e-16 * w0.norm());
        assert!((dw0 + I * FRAC_1_SQRT_PI / (z * z)).norm() <= 1e-16 * dw0.norm());
        let (w1, _) = laplace_cf(z, 1);
        let expected = I * FRAC_1_SQRT_PI * z / (z * z - 0.5);
        assert!((w1 - expected).norm() <= 1e-15 * w1.norm());
    }

    /// The derivative is that of the formula: central differences agree.
    #[test]
    fn derivative_is_that_of_the_formula() {
        for (depth, z) in [
            (20, Complex64::new(9.0, 3.0)),
            (16, Complex64::new(7.0, 0.01)),
        ] {
            let h = 1e-5 * z.norm();
            let fd = (laplace_cf(z + h, depth).0 - laplace_cf(z - h, depth).0) / (2.0 * h);
            let dw = laplace_cf(z, depth).1;
            assert!((fd - dw).norm() <= 1e-8 * dw.norm());
            let fd = (gauss_cf(z + h, depth).0 - gauss_cf(z - h, depth).0) / (2.0 * h);
            let dw = gauss_cf(z, depth).1;
            assert!((fd - dw).norm() <= 1e-8 * dw.norm());
        }
    }

    /// The derivative stays finite where 2z overflows and exp(−z²) underflows.
    #[test]
    fn gauss_cf_derivative_is_finite_for_huge_x() {
        let dw = gauss_cf(Complex64::new(1e308, 0.0), 16).1;
        assert!(dw.re.is_finite() && dw.im.is_finite(), "{dw}");
    }
}
