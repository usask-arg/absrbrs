//! Cubic Hermite interpolation of ln Q against an abscissa (T for TIPS-2021, ln T for user
//! tables), with node slopes from the second-order differences of numpy's `gradient`
//! (`edge_order=2`, valid on nonuniform grids; the secant for a two-point table). The derivative
//! returned is that of the interpolant, so Q and its derivative are exactly consistent.

/// A table of ln Q at increasing abscissae.
#[derive(Clone, Debug)]
pub(crate) struct Table {
    x: Vec<f64>,
    ln_q: Vec<f64>,
    slope: Vec<f64>,
}

impl Table {
    /// `x` strictly increasing with at least two points, `ln_q` the same length (checked by callers).
    pub(crate) fn new(x: Vec<f64>, ln_q: Vec<f64>) -> Table {
        let n = x.len();
        // Secant slopes: d[i] between nodes i and i + 1.
        let h: Vec<f64> = x.windows(2).map(|w| w[1] - w[0]).collect();
        let d: Vec<f64> = (0..n - 1).map(|i| (ln_q[i + 1] - ln_q[i]) / h[i]).collect();
        let mut slope = vec![0.0; n];
        if n == 2 {
            slope[0] = d[0];
            slope[1] = d[0];
        } else {
            // The derivatives of the quadratic through three neighbouring nodes: at the middle
            // node for the interior, at the outer node for the ends.
            for i in 1..n - 1 {
                slope[i] = (h[i - 1] * d[i] + h[i] * d[i - 1]) / (h[i - 1] + h[i]);
            }
            slope[0] = d[0] - h[0] * (d[1] - d[0]) / (h[0] + h[1]);
            slope[n - 1] = d[n - 2] + h[n - 2] * (d[n - 2] - d[n - 3]) / (h[n - 3] + h[n - 2]);
        }
        Table { x, ln_q, slope }
    }

    /// (ln Q, d ln Q/dx) at `x`, which must lie between the first and last abscissae (checked
    /// by callers).
    pub(crate) fn evaluate(&self, x: f64) -> (f64, f64) {
        let n = self.x.len();
        let i = self.x.partition_point(|&xi| xi <= x).clamp(1, n - 1) - 1;
        let h = self.x[i + 1] - self.x[i];
        let s = (x - self.x[i]) / h;
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

#[cfg(test)]
mod tests {
    use super::Table;

    #[test]
    fn slopes_are_exact_for_quadratics_on_a_nonuniform_grid() {
        let x = vec![1.0, 1.5, 3.0, 3.25, 5.0];
        let f = |x: f64| 2.0 - 0.7 * x + 0.3 * x * x;
        let table = Table::new(x.clone(), x.iter().map(|&x| f(x)).collect());
        for (xi, s) in x.iter().zip(&table.slope) {
            assert!((s - (-0.7 + 0.6 * xi)).abs() <= 1e-13, "{xi}: {s}");
        }
        let two = Table::new(vec![1.0, 3.0], vec![1.0, 2.0]);
        assert_eq!(two.slope, [0.5, 0.5]);
        assert_eq!(two.evaluate(2.0), (1.5, 0.5));
    }
}
