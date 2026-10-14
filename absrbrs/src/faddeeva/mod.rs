//! The Faddeeva function w(z) = exp(−z²) erfc(−iz) and the Voigt kernel (spec §6).
//!
//! Every function takes Im z ≥ 0 (y ≥ 0): the Voigt profile is K(x, y) = Re w(x + iy) with
//! y = γ_L/σ_D ≥ 0. Negative x uses the exact symmetries w(−x + iy) = conj(w(x + iy)) and
//! w′(−x + iy) = −conj(w′(x + iy)). For Im z < 0 or NaN the results are NaN. Accuracy is
//! validated for y ≤ 10⁶ and |x| ≤ 10⁸ (spec §6.2).

mod accurate;
mod cf;
mod fast;

pub use num_complex::Complex64;

use crate::{Error, Result};

/// 1/√π.
pub(crate) const FRAC_1_SQRT_PI: f64 = 0.564_189_583_547_756_3;

/// The accuracy tier of the Faddeeva function (spec §6.3).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum Tier {
    /// Value error ≤ 1e-6 relative to Re w; derivative error ≤ 1e-5 relative to |w′|. The
    /// region-switched composite of docs/plans/phase-0-followups.md. The default. Near the real
    /// axis (y ≤ 0.01, |x| < 6) w is a Taylor series in y, which is not analytic in z; there the
    /// derivative returned is K_x − iK_y of that formula, so its imaginary part is −K_y and not
    /// ∂ Im w/∂x.
    #[default]
    Fast,
    /// Value error ≤ 2e-13 relative to Re w; derivative error ≤ 2e-12 relative to |w′|. Value
    /// ported from the Faddeeva Package; derivative absrbrs's own. Error bounds hold over the
    /// validated domain (y ≤ 10⁶, |x| ≤ 10⁸).
    Accurate,
}

/// The Voigt kernel triplet at x + iy: K = Re w, K_x = ∂K/∂x and K_y = ∂K/∂y.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Kernel {
    /// K = Re w(x + iy).
    pub k: f64,
    /// ∂K/∂x = Re w′.
    pub k_x: f64,
    /// ∂K/∂y = −Im w′.
    pub k_y: f64,
}

/// w(z) and its derivative w′(z), for Im z ≥ 0.
///
/// For Im z < 0 or NaN, both results are NaN.
pub fn w(z: Complex64, tier: Tier) -> (Complex64, Complex64) {
    if z.im.is_nan() || z.im < 0.0 {
        let nan = Complex64::new(f64::NAN, f64::NAN);
        return (nan, nan);
    }
    let reflect = z.re < 0.0;
    // `+ 0.0` turns y = −0.0 into +0.0, so −0.0 gives the same bits as +0.0 in every region.
    let z = Complex64::new(if reflect { -z.re } else { z.re }, z.im + 0.0);
    let (w, dw) = match tier {
        Tier::Fast => fast::w_and_derivative(z),
        Tier::Accurate => accurate::w_and_derivative(z),
    };
    if reflect {
        (w.conj(), -dw.conj())
    } else {
        (w, dw)
    }
}

/// The kernel triplet at x + iy, for y ≥ 0.
///
/// For y < 0 or NaN, all three components are NaN.
pub fn kernel(x: f64, y: f64, tier: Tier) -> Kernel {
    let (w, dw) = w(Complex64::new(x, y), tier);
    Kernel {
        k: w.re,
        k_x: dw.re,
        k_y: -dw.im,
    }
}

/// The kernel triplet at every x of `x` for one y, into `k`, `k_x` and `k_y` (each as long as
/// `x`). Each result has the accuracy of [`kernel`] (the current implementation is bitwise equal
/// to it; a future batched implementation may differ in the last bits). For y < 0 or NaN, every
/// output is NaN. Returns [`Error::InvalidInput`] if the lengths differ.
pub fn kernel_into(
    x: &[f64],
    y: f64,
    tier: Tier,
    k: &mut [f64],
    k_x: &mut [f64],
    k_y: &mut [f64],
) -> Result<()> {
    let n = x.len();
    if k.len() != n || k_x.len() != n || k_y.len() != n {
        return Err(Error::InvalidInput(format!(
            "kernel_into: x has {n} values but k, k_x and k_y have {}, {} and {}",
            k.len(),
            k_x.len(),
            k_y.len()
        )));
    }
    for i in 0..n {
        let triplet = kernel(x[i], y, tier);
        (k[i], k_x[i], k_y[i]) = (triplet.k, triplet.k_x, triplet.k_y);
    }
    Ok(())
}
