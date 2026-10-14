//! Cubic Hermite interpolation of ln Q(T), with node slopes from the second-order centred
//! difference (numpy's `gradient`; first order at the ends). The derivative returned is that of
//! the interpolant, so Q and d ln Q/dT are exactly consistent.

/// A table of ln Q at increasing temperatures.
#[derive(Clone, Debug)]
pub(crate) struct Table {
    pub(crate) t: Vec<f64>,
    pub(crate) ln_q: Vec<f64>,
    slope: Vec<f64>,
}

impl Table {
    /// `t` strictly increasing with at least two points, `ln_q` the same length (checked by callers).
    pub(crate) fn new(t: Vec<f64>, ln_q: Vec<f64>) -> Table {
        let n = t.len();
        let mut slope = vec![0.0; n];
        slope[0] = (ln_q[1] - ln_q[0]) / (t[1] - t[0]);
        slope[n - 1] = (ln_q[n - 1] - ln_q[n - 2]) / (t[n - 1] - t[n - 2]);
        for i in 1..n - 1 {
            let (hl, hr) = (t[i] - t[i - 1], t[i + 1] - t[i]);
            slope[i] = (hl * hl * ln_q[i + 1] - hr * hr * ln_q[i - 1]
                + (hr * hr - hl * hl) * ln_q[i])
                / (hl * hr * (hl + hr));
        }
        Table { t, ln_q, slope }
    }

    /// The first and last temperatures.
    pub(crate) fn range(&self) -> (f64, f64) {
        (self.t[0], self.t[self.t.len() - 1])
    }

    /// (ln Q, d ln Q/dT) at `temperature`, which must lie in `range()` (checked by callers).
    pub(crate) fn evaluate(&self, temperature: f64) -> (f64, f64) {
        let n = self.t.len();
        let i = self
            .t
            .partition_point(|&t| t <= temperature)
            .clamp(1, n - 1)
            - 1;
        let h = self.t[i + 1] - self.t[i];
        let s = (temperature - self.t[i]) / h;
        let (y0, y1, m0, m1) = (
            self.ln_q[i],
            self.ln_q[i + 1],
            self.slope[i] * h,
            self.slope[i + 1] * h,
        );
        let s2 = s * s;
        let s3 = s2 * s;
        let value = (2.0 * s3 - 3.0 * s2 + 1.0) * y0
            + (s3 - 2.0 * s2 + s) * m0
            + (-2.0 * s3 + 3.0 * s2) * y1
            + (s3 - s2) * m1;
        let deriv = ((6.0 * s2 - 6.0 * s) * y0
            + (3.0 * s2 - 4.0 * s + 1.0) * m0
            + (-6.0 * s2 + 6.0 * s) * y1
            + (3.0 * s2 - 2.0 * s) * m1)
            / h;
        (value, deriv)
    }
}
