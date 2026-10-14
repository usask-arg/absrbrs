//! Weideman's (1994) N = 24 rational approximation of w(z), with the derivative of the formula.

use num_complex::Complex64;

use crate::faddeeva::FRAC_1_SQRT_PI;

/// Number of terms.
pub(crate) const N: usize = 24;
/// L = √(N/√2).
pub(crate) const L: f64 = 4.119534287814235;
/// a_N, …, a₁: highest power of Z first, for Horner's rule (reference/data/weideman.jsonl).
pub(crate) const A: [f64; N] = [
    -1.5137461285014412e-10,
    4.9048204900217325e-9,
    1.3310454420371002e-9,
    -3.008282348388502e-8,
    -1.9122258888735165e-8,
    1.8738343454069238e-7,
    2.5682641335791855e-7,
    -1.085647579569039e-6,
    -3.0388931842808516e-6,
    4.139461724333319e-6,
    3.0471066082991645e-5,
    2.4331415462260513e-5,
    -0.00020748431511424012,
    -0.0007816642995623505,
    -0.0004936426901285409,
    0.006215006362949199,
    0.033723366855315975,
    0.10838723484566723,
    0.26549639598807706,
    0.5361139535729115,
    0.9257087138588673,
    1.394819673379119,
    1.8562864992055401,
    2.1978589365315413,
];

const I: Complex64 = Complex64::new(0.0, 1.0);

/// w(z) ≈ 2 p(Z)/(L − iz)² + 1/(√π (L − iz)), Z = (L + iz)/(L − iz), p(Z) = Σ aₙ Zⁿ⁻¹, and its
/// derivative.
pub(crate) fn w_and_derivative(z: Complex64) -> (Complex64, Complex64) {
    let d = L - I * z;
    let inv = d.inv();
    let big_z = (L + I * z) * inv;
    let (mut p, mut dp) = (Complex64::new(A[0], 0.0), Complex64::new(0.0, 0.0));
    for &a in &A[1..] {
        dp = dp * big_z + p;
        p = p * big_z + a;
    }
    let inv2 = inv * inv;
    let w = 2.0 * p * inv2 + FRAC_1_SQRT_PI * inv;
    let dz = 2.0 * I * L * inv2; // dZ/dz
    let dw = 2.0 * dp * dz * inv2 + 4.0 * I * p * inv2 * inv + I * FRAC_1_SQRT_PI * inv2;
    (w, dw)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The coefficients agree with an f64 evaluation of the same discrete Fourier formula.
    #[test]
    fn coefficients_follow_from_the_formula() {
        let m = 2 * N;
        let f = |k: usize| {
            let t = L * (std::f64::consts::PI * k as f64 / (2 * m) as f64).tan();
            (-t * t).exp() * (L * L + t * t)
        };
        for n in 1..=N {
            let sum: f64 = (1..m)
                .map(|k| f(k) * (std::f64::consts::PI * (k * n) as f64 / m as f64).cos())
                .sum();
            let a = (f(0) + 2.0 * sum) / (2 * m) as f64;
            assert!((a - A[N - n]).abs() <= 1e-13, "a_{n}: {a} vs {}", A[N - n]);
        }
        assert!((L - (N as f64 / std::f64::consts::SQRT_2).sqrt()).abs() <= 1e-15 * L);
    }
}
