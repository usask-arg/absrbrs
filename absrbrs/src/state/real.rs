//! The scalar arithmetic the line-state formulas need: f64 in production, dual numbers in tests.

use std::ops::{Add, Div, Mul, Neg, Sub};

/// A real scalar: `f64`, or a dual number carrying derivatives (tests only).
#[allow(dead_code)] // used by the engine (phase 1c)
pub(crate) trait Real:
    Copy
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + Neg<Output = Self>
{
    fn from_f64(v: f64) -> Self;
    fn exp(self) -> Self;
    fn exp_m1(self) -> Self;
    fn ln(self) -> Self;
    fn sqrt(self) -> Self;
    fn powf(self, e: Self) -> Self;
}

impl Real for f64 {
    fn from_f64(v: f64) -> f64 {
        v
    }
    fn exp(self) -> f64 {
        f64::exp(self)
    }
    fn exp_m1(self) -> f64 {
        f64::exp_m1(self)
    }
    fn ln(self) -> f64 {
        f64::ln(self)
    }
    fn sqrt(self) -> f64 {
        f64::sqrt(self)
    }
    fn powf(self, e: f64) -> f64 {
        f64::powf(self, e)
    }
}
