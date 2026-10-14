//! The committed Faddeeva reference data (reference/data/faddeeva.jsonl) load, cover the grid of
//! spec §6.2 and agree with known values.

mod common;

use common::{find, load};

/// The double nearest to 2/√π.
const FRAC_2_SQRT_PI: f64 = std::f64::consts::FRAC_2_SQRT_PI;

fn close(got: f64, expected: f64) -> bool {
    (got - expected).abs() <= f64::EPSILON * expected.abs()
}

#[test]
fn faddeeva_file_covers_the_spec_grid() {
    let records = load("faddeeva.jsonl");
    assert_eq!(records.len(), 74);
    for r in &records {
        assert_eq!(r.generator, "faddeeva");
        assert_eq!(r.label, "certified");
        let x = r.input_f64s("x");
        assert_eq!(x.len(), 124);
        assert!(x.windows(2).all(|p| p[0] < p[1]));
        assert_eq!((x[0], x[x.len() - 1]), (0.0, 1e8));
        for name in ["w_re", "w_im", "dw_re", "dw_im"] {
            assert_eq!(r.output(name).hi.len(), x.len());
        }
    }
    let y: Vec<f64> = records.iter().map(|r| r.input_f64("y")).collect();
    assert!(y.windows(2).all(|p| p[0] < p[1]));
    assert_eq!((y[0], y[1], y[73]), (0.0, 1e-30, 1e6));
    assert_eq!(find(&records, "y=1.0").input_f64("y"), 1.0);
}

#[test]
fn faddeeva_known_values() {
    let records = load("faddeeva.jsonl");
    let real_axis = find(&records, "y=0.0");
    // w(0) = 1 and w'(0) = 2i/sqrt(pi), exactly representable as the nearest doubles.
    assert_eq!(real_axis.output("w_re").hi[0], 1.0);
    assert_eq!(real_axis.output("w_im").hi[0], 0.0);
    assert_eq!(real_axis.output("dw_re").hi[0], 0.0);
    assert_eq!(real_axis.output("dw_im").hi[0], FRAC_2_SQRT_PI);
    // w(1) = exp(-1) + (2i/sqrt(pi)) D(1), D Dawson's integral.
    let x = real_axis.input_f64s("x");
    let one = x
        .iter()
        .position(|&v| v == 1.0)
        .expect("x = 1 is on the grid");
    assert!(close(
        real_axis.output("w_re").hi[one],
        0.367_879_441_171_442_3
    ));
    assert!(close(
        real_axis.output("w_im").hi[one],
        0.607_157_705_841_393_7
    ));
    // w(i) = erfcx(1).
    let y1 = find(&records, "y=1.0");
    assert!(close(y1.output("w_re").hi[0], 0.427_583_576_155_807));
    assert_eq!(y1.output("w_im").hi[0], 0.0);
}
