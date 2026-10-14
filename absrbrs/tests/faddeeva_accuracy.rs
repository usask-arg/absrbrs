//! The Faddeeva tiers against the committed Arb grid (reference/data/faddeeva.jsonl), with the
//! error metrics of spec §6.2, and the API's symmetry and domain rules.

mod common;

use absrbrs::faddeeva::{self, Complex64, Kernel, Tier};
use common::load;

/// The largest value and derivative errors of `tier` on the grid, with where they occur.
fn worst(tier: Tier) -> ((f64, f64, f64), (f64, f64, f64)) {
    let (mut value, mut derivative) = ((0.0, 0.0, 0.0), (0.0, 0.0, 0.0));
    for r in load("faddeeva.jsonl") {
        let y = r.input_f64("y");
        let (w_re, dw_re, dw_im) = (r.output("w_re"), r.output("dw_re"), r.output("dw_im"));
        for (i, x) in r.input_f64s("x").into_iter().enumerate() {
            let (w, dw) = faddeeva::w(Complex64::new(x, y), tier);
            let k = w_re.hi[i];
            let ve = if k < f64::MIN_POSITIVE {
                if w.re.is_finite() { 0.0 } else { f64::INFINITY }
            } else {
                (w.re - k).abs() / k
            };
            let de = (dw.re - dw_re.hi[i]).abs().max((dw.im - dw_im.hi[i]).abs())
                / dw_re.hi[i].hypot(dw_im.hi[i]);
            let (ve, de) = (nan_is_infinite(ve), nan_is_infinite(de));
            if ve > value.0 {
                value = (ve, x, y);
            }
            if de > derivative.0 {
                derivative = (de, x, y);
            }
        }
    }
    (value, derivative)
}

fn nan_is_infinite(e: f64) -> f64 {
    if e.is_nan() { f64::INFINITY } else { e }
}

#[test]
fn accurate_tier_meets_its_targets() {
    let (v, d) = worst(Tier::Accurate);
    eprintln!(
        "Accurate: value {:e} at ({}, {}), derivative {:e} at ({}, {})",
        v.0, v.1, v.2, d.0, d.1, d.2
    );
    assert!(
        v.0 <= 2e-13,
        "value error {:e} at x = {}, y = {}",
        v.0,
        v.1,
        v.2
    );
    assert!(
        d.0 <= 2e-12,
        "derivative error {:e} at x = {}, y = {}",
        d.0,
        d.1,
        d.2
    );
}

#[test]
fn negative_x_follows_the_symmetries_exactly() {
    for tier in [Tier::Accurate] {
        for (x, y) in [
            (0.3, 0.0),
            (2.0, 0.5),
            (7.0, 1e-3),
            (40.0, 2.0),
            (1e5, 1e-20),
        ] {
            let (w, dw) = faddeeva::w(Complex64::new(x, y), tier);
            let (wn, dwn) = faddeeva::w(Complex64::new(-x, y), tier);
            assert_eq!(wn, w.conj(), "w at ({x}, {y})");
            assert_eq!(dwn, -dw.conj(), "w′ at ({x}, {y})");
        }
    }
}

#[test]
fn negative_or_nan_y_gives_nan() {
    for z in [
        Complex64::new(1.0, -1e-300),
        Complex64::new(1.0, f64::NAN),
        Complex64::new(1.0, f64::NEG_INFINITY),
    ] {
        let (w, dw) = faddeeva::w(z, Tier::Accurate);
        assert!(w.re.is_nan() && w.im.is_nan() && dw.re.is_nan() && dw.im.is_nan());
    }
}

#[test]
fn negative_zero_y_is_accepted_and_equals_positive_zero() {
    for x in [0.5, 30.0] {
        let (w, dw) = faddeeva::w(Complex64::new(x, 0.0), Tier::Accurate);
        let (wn, dwn) = faddeeva::w(Complex64::new(x, -0.0), Tier::Accurate);
        let bits = |c: Complex64| (c.re.to_bits(), c.im.to_bits());
        assert!(w.re.is_finite() && w.im.is_finite() && dw.re.is_finite() && dw.im.is_finite());
        assert_eq!(bits(wn), bits(w), "w at x = {x}");
        assert_eq!(bits(dwn), bits(dw), "w′ at x = {x}");
    }
}

#[test]
fn kernel_is_the_triplet_of_w() {
    for (x, y) in [(0.0, 0.0), (1.5, 0.25), (-9.0, 1e-3), (300.0, 4.0)] {
        let (w, dw) = faddeeva::w(Complex64::new(x, y), Tier::Accurate);
        assert_eq!(
            faddeeva::kernel(x, y, Tier::Accurate),
            Kernel {
                k: w.re,
                k_x: dw.re,
                k_y: -dw.im
            }
        );
    }
}

#[test]
fn kernel_into_matches_kernel_bitwise() {
    let x: Vec<f64> = (-400..=400).map(|i| f64::from(i) * 0.137).collect();
    for y in [0.0, 1e-3, 0.5, 30.0] {
        let n = x.len();
        let (mut k, mut kx, mut ky) = (vec![0.0; n], vec![0.0; n], vec![0.0; n]);
        faddeeva::kernel_into(&x, y, Tier::Accurate, &mut k, &mut kx, &mut ky).unwrap();
        for i in 0..n {
            let one = faddeeva::kernel(x[i], y, Tier::Accurate);
            assert_eq!(
                (k[i], kx[i], ky[i]),
                (one.k, one.k_x, one.k_y),
                "x = {}, y = {y}",
                x[i]
            );
        }
    }
}

#[test]
fn kernel_into_rejects_mismatched_lengths() {
    let (mut k, mut kx, mut ky) = (vec![0.0; 3], vec![0.0; 3], vec![0.0; 2]);
    let err = faddeeva::kernel_into(
        &[0.0, 1.0, 2.0],
        1.0,
        Tier::Accurate,
        &mut k,
        &mut kx,
        &mut ky,
    );
    assert!(matches!(err, Err(absrbrs::Error::InvalidInput(_))));
}
