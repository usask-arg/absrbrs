//! The Faddeeva function w(z) = exp(−z²) erfc(−iz), for complex z.
//!
//! Ported from `Faddeeva.cc` (revision of 12 May 2015): `sinc`, `sinh_taylor`, `sqr`, `expa2n2`
//! and `FADDEEVA(w)`. Only the machine-precision path (`relerr` = `DBL_EPSILON`) and the
//! `USE_CONTINUED_FRACTION` branch are ported.
// Copyright (c) 2012 Massachusetts Institute of Technology
//
// Permission is hereby granted, free of charge, to any person obtaining
// a copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#![allow(clippy::excessive_precision)] // the coefficients are copied verbatim

// To compute the Faddeeva function, we use a combination of two
// algorithms:
//
// For sufficiently large |z|, we use a continued-fraction expansion
// for w(z) similar to those described in:
//
//    Walter Gautschi, "Efficient computation of the complex error
//    function," SIAM J. Numer. Anal. 7(1), pp. 187-198 (1970)
//
//    G. P. M. Poppe and C. M. J. Wijers, "More efficient computation
//    of the complex error function," ACM Trans. Math. Soft. 16(1),
//    pp. 38-46 (1990).
//
// Unlike those papers, however, we switch to a completely different
// algorithm for smaller |z|:
//
//    Mofreh R. Zaghloul and Ahmed N. Ali, "Algorithm 916: Computing the
//    Faddeyeva and Voigt Functions," ACM Trans. Math. Soft. 38(2), 15
//    (2011).
//
// (I initially used this algorithm for all z, but it turned out to be
//  significantly slower than the continued-fraction expansion for
//  larger |z|.  On the other hand, it is competitive for smaller |z|,
//  and is significantly more accurate than the Poppe & Wijers code
//  in some regions, e.g. in the vicinity of z=1+1i.)
//
// Note that this is an INDEPENDENT RE-IMPLEMENTATION of these algorithms,
// based on the description in the papers ONLY.  In particular, I did
// not refer to the authors' Fortran or Matlab implementations, respectively,
// (which are under restrictive ACM copyright terms and therefore unusable
//  in free/open-source software).
//
// Steven G. Johnson, Massachusetts Institute of Technology
// http://math.mit.edu/~stevenj
// October 2012.
//
//  -- Note that Algorithm 916 assumes that the erfc(x) function,
//     or rather the scaled function erfcx(x) = exp(x*x)*erfc(x),
//     is supplied for REAL arguments x.   I originally used an
//     erfcx routine derived from DERFC in SLATEC, but I have
//     since replaced it with a much faster routine written by
//     me which uses a combination of continued-fraction expansions
//     and a lookup table of Chebyshev polynomials.  For speed,
//     I implemented a similar algorithm for Im[w(x)] of real x,
//     since this comes up frequently in the other error functions.

use num_complex::Complex64;

// return sinc(x) = sin(x)/x, given both x and sin(x)
// [since we only use this in cases where sin(x) has already been computed]
#[inline]
fn sinc(x: f64, sinx: f64) -> f64 {
    if x.abs() < 1e-4 {
        1.0 - (0.1666666666666666666667) * x * x
    } else {
        sinx / x
    }
}

// sinh(x) via Taylor series, accurate to machine precision for |x| < 1e-2
#[inline]
fn sinh_taylor(x: f64) -> f64 {
    x * (1.0 + (x * x) * (0.1666666666666666666667 + 0.00833333333333333333333 * (x * x)))
}

#[inline]
fn sqr(x: f64) -> f64 {
    x * x
}

// precomputed table of expa2n2[n-1] = exp(-a2*n*n)
// for double-precision a2 = 0.26865... in FADDEEVA(w), below.
#[rustfmt::skip] // one line per line of the C table, for comparison with Faddeeva.cc
const EXPA2N2: [f64; 52] = [
    7.64405281671221563e-01,
    3.41424527166548425e-01,
    8.91072646929412548e-02,
    1.35887299055460086e-02,
    1.21085455253437481e-03,
    6.30452613933449404e-05,
    1.91805156577114683e-06,
    3.40969447714832381e-08,
    3.54175089099469393e-10,
    2.14965079583260682e-12,
    7.62368911833724354e-15,
    1.57982797110681093e-17,
    1.91294189103582677e-20,
    1.35344656764205340e-23,
    5.59535712428588720e-27,
    1.35164257972401769e-30,
    1.90784582843501167e-34,
    1.57351920291442930e-38,
    7.58312432328032845e-43,
    2.13536275438697082e-47,
    3.51352063787195769e-52,
    3.37800830266396920e-57,
    1.89769439468301000e-62,
    6.22929926072668851e-68,
    1.19481172006938722e-73,
    1.33908181133005953e-79,
    8.76924303483223939e-86,
    3.35555576166254986e-92,
    7.50264110688173024e-99,
    9.80192200745410268e-106,
    7.48265412822268959e-113,
    3.33770122566809425e-120,
    8.69934598159861140e-128,
    1.32486951484088852e-135,
    1.17898144201315253e-143,
    6.13039120236180012e-152,
    1.86258785950822098e-160,
    3.30668408201432783e-169,
    3.43017280887946235e-178,
    2.07915397775808219e-187,
    7.36384545323984966e-197,
    1.52394760394085741e-206,
    1.84281935046532100e-216,
    1.30209553802992923e-226,
    5.37588903521080531e-237,
    1.29689584599763145e-247,
    1.82813078022866562e-258,
    1.50576355348684241e-269,
    7.24692320799294194e-281,
    2.03797051314726829e-292,
    3.34880215927873807e-304,
    0.0, // underflow (also prevents reads past array end, below)
];

/// w(z) to machine precision (Johnson's `relerr` = `DBL_EPSILON`).
///
/// The y < 0 branches are ported for comparison with the source; callers in this crate pass
/// Re z ≥ 0 and Im z ≥ 0 only.
pub(crate) fn w(z: Complex64) -> Complex64 {
    if z.re == 0.0 {
        return Complex64::new(super::erfcx(z.im), z.re); // give correct sign of 0 in cimag(w)
    } else if z.im == 0.0 {
        return Complex64::new(f64::exp(-sqr(z.re)), super::w_im(z.re));
    }

    // The source takes `relerr` as a parameter and computes a, c and a2 from it when
    // relerr > DBL_EPSILON. Only the `relerr <= DBL_EPSILON` branch is ported: relerr is
    // DBL_EPSILON and a, c, a2 are the precomputed constants.
    let relerr = f64::EPSILON;
    let a: f64 = 0.518321480430085929872; // pi / sqrt(-log(eps*0.5))
    let c: f64 = 0.329973702884629072537; // (2/pi) * a;
    let a2: f64 = 0.268657157075235951582; // a^2
    let x = z.re.abs();
    let y = z.im;
    let ya = y.abs();

    // return value (`cmplx ret = 0.` in the source; every path below assigns it)
    let ret: Complex64;

    let (mut sum1, mut sum2, mut sum3, mut sum4, mut sum5) = (0.0, 0.0, 0.0, 0.0, 0.0);

    // USE_CONTINUED_FRACTION is 1 in the source: only that branch is ported.
    if ya > 7.0
        || (x > 6.0 // continued fraction is faster
            // As pointed out by M. Zaghloul, the continued
            // fraction seems to give a large relative error in
            // Re w(z) for |x| ~ 6 and small |y|, so use
            // algorithm 816 in this region:
            && (ya > 0.1 || (x > 8.0 && ya > 1e-10) || x > 28.0))
    {
        // Poppe & Wijers suggest using a number of terms
        //     nu = 3 + 1442 / (26*rho + 77)
        // where rho = sqrt((x/x0)^2 + (y/y0)^2) where x0=6.3, y0=4.4.
        // (They only use this expansion for rho >= 1, but rho a little less
        //  than 1 seems okay too.)
        // Instead, I did my own fit to a slightly different function
        // that avoids the hypotenuse calculation, using NLopt to minimize
        // the sum of the squares of the errors in nu with the constraint
        // that the estimated nu be >= minimum nu to attain machine precision.
        // I also separate the regions where nu == 2 and nu == 1.
        const ISPI: f64 = 0.56418958354775628694807945156; // 1 / sqrt(pi)
        let xs = if y < 0.0 { -z.re } else { z.re }; // compute for -z if y < 0
        if x + ya > 4000.0 {
            // nu <= 2
            if x + ya > 1e7 {
                // nu == 1, w(z) = i/sqrt(pi) / z
                // scale to avoid overflow
                if x > ya {
                    let yax = ya / xs;
                    let denom = ISPI / (xs + yax * ya);
                    ret = Complex64::new(denom * yax, denom);
                } else if ya.is_infinite() {
                    return if x.is_nan() || y < 0.0 {
                        Complex64::new(f64::NAN, f64::NAN)
                    } else {
                        Complex64::new(0.0, 0.0)
                    };
                } else {
                    let xya = xs / ya;
                    let denom = ISPI / (xya * xs + ya);
                    ret = Complex64::new(denom, denom * xya);
                }
            } else {
                // nu == 2, w(z) = i/sqrt(pi) * z / (z*z - 0.5)
                let dr = xs * xs - ya * ya - 0.5;
                let di = 2.0 * xs * ya;
                let denom = ISPI / (dr * dr + di * di);
                ret = Complex64::new(denom * (xs * di - ya * dr), denom * (xs * dr + ya * di));
            }
        } else {
            // compute nu(z) estimate and do general continued fraction
            let (c0, c1, c2, c3, c4) = (3.9, 11.398, 0.08254, 0.1421, 0.2023); // fit
            let mut nu = f64::floor(c0 + c1 / (c2 * x + c3 * ya + c4));
            let mut wr = xs;
            let mut wi = ya;
            nu = 0.5 * (nu - 1.0);
            while nu > 0.4 {
                // w <- z - nu/w:
                let denom = nu / (wr * wr + wi * wi);
                wr = xs - wr * denom;
                wi = ya + wi * denom;
                nu -= 0.5;
            }
            {
                // w(z) = i/sqrt(pi) / w:
                let denom = ISPI / (wr * wr + wi * wi);
                ret = Complex64::new(denom * wi, denom * wr);
            }
        }
        if y < 0.0 {
            // use w(z) = 2.0*exp(-z*z) - w(-z),
            // but be careful of overflow in exp(-z*z)
            //                                = exp(-(xs*xs-ya*ya) -2*i*xs*ya)
            return 2.0 * Complex64::new((ya - xs) * (xs + ya), 2.0 * xs * y).exp() - ret;
        } else {
            return ret;
        }
    } else if x < 10.0 {
        // Note: The test that seems to be suggested in the paper is x <
        // sqrt(-log(DBL_MIN)), about 26.6, since otherwise exp(-x^2)
        // underflows to zero and sum1,sum2,sum4 are zero.  However, long
        // before this occurs, the sum1,sum2,sum4 contributions are
        // negligible in double precision; I find that this happens for x >
        // about 6, for all y.  On the other hand, I find that the case
        // where we compute all of the sums is faster (at least with the
        // precomputed expa2n2 table) until about x=10.  Furthermore, if we
        // try to compute all of the sums for x > 20, I find that we
        // sometimes run into numerical problems because underflow/overflow
        // problems start to appear in the various coefficients of the sums,
        // below.  Therefore, we use x < 10 here.
        let mut prod2ax = 1.0;
        let mut prodm2ax = 1.0;
        let expx2;

        if y.is_nan() {
            return Complex64::new(y, y);
        }

        // Somewhat ugly copy-and-paste duplication here, but I see significant
        // speedups from using the special-case code with the precomputed
        // exponential, and the x < 5e-4 special case is needed for accuracy.

        // The source tests `relerr == DBL_EPSILON` here, which always holds in this port: only
        // the branch that uses the precomputed exp(-a2*(n*n)) table is ported.
        if x < 5e-4 {
            // compute sum4 and sum5 together as sum5-sum4
            let x2 = x * x;
            expx2 = 1.0 - x2 * (1.0 - 0.5 * x2); // exp(-x*x) via Taylor
            // compute exp(2*a*x) and exp(-2*a*x) via Taylor, to double precision
            let ax2 = 1.036642960860171859744 * x; // 2*a*x
            let exp2ax = 1.0 + ax2 * (1.0 + ax2 * (0.5 + 0.166666666666666666667 * ax2));
            let expm2ax = 1.0 - ax2 * (1.0 - ax2 * (0.5 - 0.166666666666666666667 * ax2));
            let mut n: i32 = 1;
            loop {
                // The table index is in range: EXPA2N2[51] is 0.0, so coef is 0 there and the
                // convergence test below stops the loop (sum3 > 0).
                let coef = EXPA2N2[(n - 1) as usize] * expx2 / (a2 * f64::from(n * n) + y * y);
                prod2ax *= exp2ax;
                prodm2ax *= expm2ax;
                sum1 += coef;
                sum2 += coef * prodm2ax;
                sum3 += coef * prod2ax;

                // really = sum5 - sum4
                sum5 += coef * (2.0 * a) * f64::from(n) * sinh_taylor((2.0 * a) * f64::from(n) * x);

                // test convergence via sum3
                if coef * prod2ax < relerr * sum3 {
                    break;
                }
                n += 1;
            }
        } else {
            // x > 5e-4, compute sum4 and sum5 separately
            expx2 = f64::exp(-x * x);
            let exp2ax = f64::exp((2.0 * a) * x);
            let expm2ax = 1.0 / exp2ax;
            let mut n: i32 = 1;
            loop {
                // In range as above (sum5 > 0).
                let coef = EXPA2N2[(n - 1) as usize] * expx2 / (a2 * f64::from(n * n) + y * y);
                prod2ax *= exp2ax;
                prodm2ax *= expm2ax;
                sum1 += coef;
                sum2 += coef * prodm2ax;
                sum4 += (coef * prodm2ax) * (a * f64::from(n));
                sum3 += coef * prod2ax;
                sum5 += (coef * prod2ax) * (a * f64::from(n));
                // test convergence via sum5, since this sum has the slowest decay
                if (coef * prod2ax) * (a * f64::from(n)) < relerr * sum5 {
                    break;
                }
                n += 1;
            }
        }
        let expx2erfcxy = // avoid spurious overflow for large negative y
            if y > -6.0 { // for y < -6, erfcx(y) = 2*exp(y*y) to double precision
                expx2 * super::erfcx(y)
            } else {
                2.0 * f64::exp(y * y - x * x)
            };
        if y > 5.0 {
            // imaginary terms cancel
            let sinxy = f64::sin(x * y);
            ret = Complex64::new(
                (expx2erfcxy - c * y * sum1) * f64::cos(2.0 * x * y)
                    + (c * x * expx2) * sinxy * sinc(x * y, sinxy),
                0.0,
            );
        } else {
            let xs = z.re;
            let sinxy = f64::sin(xs * y);
            let sin2xy = f64::sin(2.0 * xs * y);
            let cos2xy = f64::cos(2.0 * xs * y);
            let coef1 = expx2erfcxy - c * y * sum1;
            let coef2 = c * xs * expx2;
            ret = Complex64::new(
                coef1 * cos2xy + coef2 * sinxy * sinc(xs * y, sinxy),
                coef2 * sinc(2.0 * xs * y, sin2xy) - coef1 * sin2xy,
            );
        }
    } else {
        // x large: only sum3 & sum5 contribute (see above note)
        if x.is_nan() {
            return Complex64::new(x, x);
        }
        if y.is_nan() {
            return Complex64::new(y, y);
        }

        // (USE_CONTINUED_FRACTION branch)
        ret = Complex64::new(f64::exp(-x * x), 0.0); // |y| < 1e-10, so we only need exp(-x*x) term
        // (round instead of ceil as in original paper; note that x/a > 1 here)
        let n0 = f64::floor(x / a + 0.5); // sum in both directions, starting at n0
        let dx = a * n0 - x;
        sum3 = f64::exp(-dx * dx) / (a2 * (n0 * n0) + y * y);
        sum5 = a * n0 * sum3;
        let exp1 = f64::exp(4.0 * a * dx);
        let mut exp1dn = 1.0;
        let mut dn: i32 = 1;
        // The labelled block stands for the source's `goto finish`.
        'finish: {
            while n0 - f64::from(dn) > 0.0 {
                // loop over n0-dn and n0+dn terms
                let np = n0 + f64::from(dn);
                let nm = n0 - f64::from(dn);
                let mut tp = f64::exp(-sqr(a * f64::from(dn) + dx));
                exp1dn *= exp1;
                let mut tm = tp * exp1dn; // trick to get tm from tp
                tp /= a2 * (np * np) + y * y;
                tm /= a2 * (nm * nm) + y * y;
                sum3 += tp + tm;
                sum5 += a * (np * tp + nm * tm);
                if a * (np * tp + nm * tm) < relerr * sum5 {
                    break 'finish;
                }
                dn += 1;
            }
            loop {
                // loop over n0+dn terms only (since n0-dn <= 0)
                // `n0 + dn++` in the source: np uses dn before the increment and the exponent
                // after it, as there.
                let np = n0 + f64::from(dn);
                dn += 1;
                let tp = f64::exp(-sqr(a * f64::from(dn) + dx)) / (a2 * (np * np) + y * y);
                sum3 += tp;
                sum5 += a * np * tp;
                if a * np * tp < relerr * sum5 {
                    break 'finish;
                }
            }
        }
    }
    // finish:
    ret + Complex64::new(
        (0.5 * c) * y * (sum2 + sum3),
        (0.5 * c) * f64::copysign(sum5 - sum4, z.re),
    )
}
