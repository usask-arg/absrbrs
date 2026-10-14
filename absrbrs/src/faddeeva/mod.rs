//! The Faddeeva function w(z) = exp(−z²) erfc(−iz) and the Voigt kernel (spec §6).

#[allow(dead_code, unused_imports)] // used by the public API from Task 2 on
mod accurate;

/// 1/√π.
#[allow(dead_code)] // used from Task 2 on
pub(crate) const FRAC_1_SQRT_PI: f64 = 0.564_189_583_547_756_3;
