//! Forward-mode dual numbers with `N` derivative components, for checking the closed forms.

use std::ops::{Add, Div, Mul, Neg, Sub};

use super::Real;

/// v + Σₖ d[k] εₖ, with εⱼ εₖ = 0.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Dual<const N: usize> {
    pub v: f64,
    pub d: [f64; N],
}

impl<const N: usize> Dual<N> {
    /// A value with no dependence on the variables.
    pub fn constant(v: f64) -> Self {
        Dual { v, d: [0.0; N] }
    }

    /// The variable number `k`: d[k] = 1.
    pub fn variable(v: f64, k: usize) -> Self {
        let mut d = [0.0; N];
        d[k] = 1.0;
        Dual { v, d }
    }

    /// (f(v), f′(v) d): the chain rule for a function of one argument.
    fn chain(self, v: f64, dv: f64) -> Self {
        Dual {
            v,
            d: self.d.map(|d| dv * d),
        }
    }
}

impl<const N: usize> Add for Dual<N> {
    type Output = Self;
    fn add(self, o: Self) -> Self {
        Dual {
            v: self.v + o.v,
            d: std::array::from_fn(|k| self.d[k] + o.d[k]),
        }
    }
}

impl<const N: usize> Sub for Dual<N> {
    type Output = Self;
    fn sub(self, o: Self) -> Self {
        Dual {
            v: self.v - o.v,
            d: std::array::from_fn(|k| self.d[k] - o.d[k]),
        }
    }
}

#[allow(clippy::suspicious_arithmetic_impl)] // the product rule
impl<const N: usize> Mul for Dual<N> {
    type Output = Self;
    fn mul(self, o: Self) -> Self {
        Dual {
            v: self.v * o.v,
            d: std::array::from_fn(|k| self.d[k] * o.v + self.v * o.d[k]),
        }
    }
}

#[allow(clippy::suspicious_arithmetic_impl)] // the quotient rule
impl<const N: usize> Div for Dual<N> {
    type Output = Self;
    fn div(self, o: Self) -> Self {
        let v = self.v / o.v;
        Dual {
            v,
            d: std::array::from_fn(|k| (self.d[k] - v * o.d[k]) / o.v),
        }
    }
}

impl<const N: usize> Neg for Dual<N> {
    type Output = Self;
    fn neg(self) -> Self {
        Dual {
            v: -self.v,
            d: self.d.map(|d| -d),
        }
    }
}

impl<const N: usize> Real for Dual<N> {
    fn from_f64(v: f64) -> Self {
        Dual::constant(v)
    }
    fn exp(self) -> Self {
        let e = self.v.exp();
        self.chain(e, e)
    }
    fn exp_m1(self) -> Self {
        self.chain(self.v.exp_m1(), self.v.exp())
    }
    fn ln(self) -> Self {
        self.chain(self.v.ln(), 1.0 / self.v)
    }
    fn sqrt(self) -> Self {
        let s = self.v.sqrt();
        self.chain(s, 1.0 / (2.0 * s))
    }
    fn powf(self, e: Self) -> Self {
        let w = self.v.powf(e.v);
        let ln_v = self.v.ln();
        Dual {
            v: w,
            d: std::array::from_fn(|k| w * (e.d[k] * ln_v + e.v * self.d[k] / self.v)),
        }
    }
}
