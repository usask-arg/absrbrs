//! The Accurate tier: a Rust port of Steven G. Johnson's Faddeeva Package (`Faddeeva.cc`,
//! revision of 12 May 2015), MIT licence; see `absrbrs/LICENSES/Faddeeva-MIT.txt`. The
//! derivative is absrbrs's own (`w_and_derivative`).

mod erfcx;
mod w;
mod w_im;

use num_complex::Complex64;

pub(crate) use erfcx::erfcx;
pub(crate) use w_im::w_im;

use super::FRAC_1_SQRT_PI;
use super::cf::{gauss_cf, laplace_cf};

/// w(z) and w′(z) for Re z ≥ 0 and Im z ≥ 0.
///
/// The value is Johnson's. The derivative uses the identity w′ = −2zw + 2i/√π for |z| < 6,
/// where its cancellation (about |z|² in relative terms) stays within the tier's target. Beyond
/// that it uses the derivative of an asymptotic formula: exp(−z²) plus a depth-16 continued
/// fraction near the real axis (x ≥ 6, y ≤ 0.1), and a depth-20 continued fraction elsewhere.
/// Measured against Arb before implementation: 9.3e-13, 2.8e-13 and 1.7e-13 relative to |w′|
/// in the three regions.
pub(crate) fn w_and_derivative(z: Complex64) -> (Complex64, Complex64) {
    let w = w::w(z);
    let dw = if z.re >= 6.0 && z.im <= 0.1 {
        gauss_cf(z, 16).1
    } else if z.norm() >= 6.0 {
        laplace_cf(z, 20).1
    } else {
        -2.0 * z * w + Complex64::new(0.0, 2.0 * FRAC_1_SQRT_PI)
    };
    (w, dw)
}
