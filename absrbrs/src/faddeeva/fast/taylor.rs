//! The small-y Taylor series of w about the real axis, with the formula's own derivatives.

use num_complex::Complex64;

use crate::faddeeva::{FRAC_1_SQRT_PI, accurate};

/// Order of the series.
const ORDER: usize = 3;

/// w(x + iy) ≈ Σₖ (iy)ᵏ w⁽ᵏ⁾(x)/k!, k = 0 … ORDER, with w(x) = exp(−x²) + i w_im(x) and
/// w⁽ᵏ⁺¹⁾ = −2x w⁽ᵏ⁾ − 2k w⁽ᵏ⁻¹⁾. The series is not analytic in z, so the derivative returned is
/// K_x − i K_y of the formula: K_x = Re Σₖ (iy)ᵏ w⁽ᵏ⁺¹⁾/k! and K_y = Re(i Σₖ₍ₖ<ORDER₎ (iy)ᵏ w⁽ᵏ⁺¹⁾/k!).
/// For |x| < 6 only: w′(x) comes from the identity −2xw + 2i/√π, which amplifies the error of
/// w_im by about 2x² (spec §6.1; accepted because w_im is accurate to about 1e-15).
pub(crate) fn w_and_derivative(x: f64, y: f64) -> (Complex64, Complex64) {
    let mut d = [Complex64::new(0.0, 0.0); ORDER + 2];
    d[0] = Complex64::new((-x * x).exp(), accurate::w_im(x));
    d[1] = -2.0 * x * d[0] + Complex64::new(0.0, 2.0 * FRAC_1_SQRT_PI);
    for k in 1..=ORDER {
        d[k + 1] = -2.0 * x * d[k] - 2.0 * k as f64 * d[k - 1];
    }
    let iy = Complex64::new(0.0, y);
    let zero = Complex64::new(0.0, 0.0);
    let (mut w, mut d_x, mut d_y) = (zero, zero, zero);
    let mut c = Complex64::new(1.0, 0.0); // (iy)ᵏ/k!
    for k in 0..=ORDER {
        w += c * d[k];
        d_x += c * d[k + 1];
        if k < ORDER {
            d_y += c * d[k + 1];
        }
        c = c * iy / (k + 1) as f64;
    }
    let k_x = d_x.re;
    let k_y = (Complex64::new(0.0, 1.0) * d_y).re;
    (w, Complex64::new(k_x, -k_y))
}
