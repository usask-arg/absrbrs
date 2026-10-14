//! The Fast tier against the Accurate tier (itself within 2e-13 of Arb): on random points, and
//! on both sides of every region boundary of the Fast composite.

use absrbrs::faddeeva::{self, Complex64, Tier};
use proptest::prelude::*;

/// Value and derivative errors of the Fast tier at z, with the Accurate tier as reference.
fn errors(z: Complex64) -> (f64, f64) {
    let (w, dw) = faddeeva::w(z, Tier::Fast);
    let (wr, dwr) = faddeeva::w(z, Tier::Accurate);
    let ve = if wr.re < f64::MIN_POSITIVE {
        if w.re.is_finite() { 0.0 } else { f64::INFINITY }
    } else {
        (w.re - wr.re).abs() / wr.re
    };
    let de = (dw.re - dwr.re).abs().max((dw.im - dwr.im).abs()) / dwr.norm();
    (
        if ve.is_nan() { f64::INFINITY } else { ve },
        if de.is_nan() { f64::INFINITY } else { de },
    )
}

fn x_values() -> impl Strategy<Value = f64> {
    prop_oneof![
        0.0..30.0f64,
        0.0..12.0f64,
        (-8.0..8.0f64).prop_map(|e| 10f64.powf(e))
    ]
}

fn y_values() -> impl Strategy<Value = f64> {
    prop_oneof![
        1 => Just(0.0),
        3 => (-4.0..0.0f64).prop_map(|e| 10f64.powf(e)),
        6 => (-30.0..6.0f64).prop_map(|e| 10f64.powf(e)),
    ]
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(20_000))]

    #[test]
    fn fast_tier_is_within_its_targets(x in x_values(), y in y_values(), negative in any::<bool>()) {
        let x = if negative { -x } else { x };
        let (ve, de) = errors(Complex64::new(x, y));
        prop_assert!(ve <= 1e-6, "value error {ve:e} at ({x}, {y})");
        prop_assert!(de <= 1e-5, "derivative error {de:e} at ({x}, {y})");
    }
}

/// The next double above a positive v.
fn up(v: f64) -> f64 {
    f64::from_bits(v.to_bits() + 1)
}

/// The next double below a positive v.
fn down(v: f64) -> f64 {
    f64::from_bits(v.to_bits() - 1)
}

/// Points on both sides of every boundary of the composite: the |z| circles at all angles and
/// their near-axis corners, and the x and y thresholds of the near-axis regions.
#[test]
fn fast_tier_is_within_its_targets_at_every_boundary() {
    let mut points = Vec::new();
    for r in [8.0, 9.0, 15.0, 50.0, 1750.0f64] {
        for j in 0..=200 {
            let theta = std::f64::consts::FRAC_PI_2 * f64::from(j) / 200.0;
            for s in [1.0 - 1e-15, 1.0, 1.0 + 1e-15] {
                points.push((r * s * theta.cos(), r * s * theta.sin()));
            }
        }
    }
    for r in [8.0, 9.0, 15.0, 50.0, 1750.0f64] {
        for y in [1e-300, 1e-30, 1e-10, 1e-4, 0.01, up(0.01)] {
            points.extend([(down(r), y), (r, y)]);
        }
    }
    for x in [6.0, 27.0f64] {
        for y in [0.0, 1e-30, 1e-10, 1e-4, 0.01, up(0.01)] {
            points.extend([(down(x), y), (x, y)]);
        }
    }
    for i in 0..=600 {
        let x = f64::from(i) * 0.05;
        points.extend([(x, 0.01), (x, up(0.01))]);
    }
    let mut worst = (0.0f64, 0.0f64);
    for (x, y) in points {
        let (ve, de) = errors(Complex64::new(x, y));
        assert!(ve <= 1e-6, "value error {ve:e} at ({x}, {y})");
        assert!(de <= 1e-5, "derivative error {de:e} at ({x}, {y})");
        worst = (worst.0.max(ve), worst.1.max(de));
    }
    eprintln!(
        "boundaries: worst value {:e}, derivative {:e}",
        worst.0, worst.1
    );
}
