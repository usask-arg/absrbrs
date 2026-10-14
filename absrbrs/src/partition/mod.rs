//! Partition functions Q(T) and d ln Q/dT (spec §5.5).

mod hermite;
mod tabulated;
mod tips2021;

pub use tabulated::TabulatedPartition;
pub use tips2021::{Tips2021, tips2021};

use crate::{Error, Result};

/// Q(T) and its logarithmic derivative d ln Q/dT (K⁻¹), consistent with each other.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PartitionValue {
    /// The total internal partition sum Q(T).
    pub q: f64,
    /// d ln Q/dT, K⁻¹.
    pub dlnq_dt: f64,
}

/// A source of partition functions, by HITRAN molecule and isotopologue number.
pub trait PartitionFunction: Send + Sync {
    /// Q and d ln Q/dT at `temperature` (K). Errors: `MissingPartitionFunction` if there is no
    /// table for the isotopologue, `OutsideTable` if the temperature is outside it (or NaN).
    fn evaluate(
        &self,
        molecule: u16,
        isotopologue: u16,
        temperature: f64,
    ) -> Result<PartitionValue>;
}

/// Returns `OutsideTable` unless `temperature` lies in `[t_first, t_last]` (so NaN is outside).
fn check_range(
    (t_first, t_last): (f64, f64),
    molecule: u16,
    isotopologue: u16,
    temperature: f64,
) -> Result<()> {
    if t_first <= temperature && temperature <= t_last {
        Ok(())
    } else {
        Err(Error::OutsideTable {
            molecule,
            isotopologue,
            temperature,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Error;

    #[test]
    fn tips_reproduces_the_source_at_nodes() {
        // Values from the Zenodo .QTpy files (1_1: Q(296) = 174.58130; 2_1: 286.09382;
        // 7_1: 215.734504; 3_1: 3474.9978), within f32 rounding.
        for (mol, iso, q296) in [
            (1, 1, 174.58130),
            (2, 1, 286.09382),
            (7, 1, 215.734504),
            (3, 1, 3474.9978),
        ] {
            let v = tips2021().evaluate(mol, iso, 296.0).unwrap();
            assert!(
                (v.q - q296).abs() <= 6e-8 * q296,
                "({mol}, {iso}): {} vs {q296}",
                v.q
            );
        }
        assert_eq!(tips2021().evaluate(1, 1, 1.0).unwrap().q, 1.0);
    }

    #[test]
    fn tips_derivative_is_that_of_the_interpolant() {
        for (mol, iso) in [(1, 1), (2, 1), (7, 1), (3, 1), (2, 10)] {
            for t in [3.7, 19.5, 150.25, 296.0, 299.9, 777.7] {
                let h = 1e-4;
                let lp = tips2021().evaluate(mol, iso, t + h).unwrap().q.ln();
                let lm = tips2021().evaluate(mol, iso, t - h).unwrap().q.ln();
                let d = tips2021().evaluate(mol, iso, t).unwrap().dlnq_dt;
                assert!(
                    ((lp - lm) / (2.0 * h) - d).abs() <= 1e-6 * d.abs().max(1e-3),
                    "({mol}, {iso}) at {t}"
                );
            }
        }
    }

    #[test]
    fn tips_coverage_and_errors() {
        let isos = tips2021().isotopologues();
        assert_eq!(isos.len(), 181);
        assert_eq!(tips2021().temperature_range(1, 1), Some((1.0, 5000.0)));
        assert_eq!(tips2021().temperature_range(3, 1), Some((1.0, 1000.0)));
        assert!(matches!(
            tips2021().evaluate(3, 1, 1000.5),
            Err(Error::OutsideTable {
                molecule: 3,
                isotopologue: 1,
                ..
            })
        ));
        assert!(matches!(
            tips2021().evaluate(1, 1, 0.5),
            Err(Error::OutsideTable { .. })
        ));
        assert!(matches!(
            tips2021().evaluate(1, 1, f64::NAN),
            Err(Error::OutsideTable { .. })
        ));
        assert!(matches!(
            tips2021().evaluate(10, 3, 296.0),
            Err(Error::MissingPartitionFunction {
                molecule: 10,
                isotopologue: 3
            })
        ));
        // Every table decodes to finite, positive Q. The converter starts the H2S tables 31_2 and
        // 31_3 at 30 K (unphysical below) and 57_1 at 2 K (Q(1 K) = 0).
        assert_eq!(tips2021().temperature_range(31, 2), Some((30.0, 5000.0)));
        assert_eq!(tips2021().temperature_range(31, 3), Some((30.0, 5000.0)));
        assert_eq!(tips2021().temperature_range(57, 1), Some((2.0, 5000.0)));
        for (mol, iso) in isos {
            let (t_first, t_last) = tips2021().temperature_range(mol, iso).unwrap();
            for t in [t_first, t_first + 0.5, 296.0f64.min(t_last), t_last] {
                let v = tips2021().evaluate(mol, iso, t).unwrap();
                assert!(
                    v.q.is_finite() && v.q > 0.0 && v.dlnq_dt.is_finite(),
                    "({mol}, {iso}) at {t}: {v:?}"
                );
            }
        }
    }

    #[test]
    fn tabulated_tables() {
        let mut p = TabulatedPartition::new();
        // Q = T^1.5 on a nonuniform grid: ln Q is linear in ln T, which the interpolant
        // reproduces to rounding, at the nodes and between them.
        let t: Vec<f64> = vec![100.0, 150.0, 220.0, 300.0, 400.0];
        let q: Vec<f64> = t.iter().map(|t: &f64| t.powf(1.5)).collect();
        p.insert(2, 1, &t, &q).unwrap();
        // A two-point table (secant slopes) is exact for a power law too.
        p.insert(2, 4, &[100.0, 400.0], &[1000.0, 8000.0]).unwrap();
        for (iso, temperature) in [
            (1, 100.0),
            (1, 123.4),
            (1, 220.0),
            (1, 250.0),
            (1, 399.9),
            (1, 400.0),
            (4, 250.0),
        ] {
            let v = p.evaluate(2, iso, temperature).unwrap();
            assert!(
                (v.q / f64::powf(temperature, 1.5) - 1.0).abs() <= 1e-12,
                "{iso} at {temperature}: {v:?}"
            );
            assert!(
                (v.dlnq_dt * temperature / 1.5 - 1.0).abs() <= 1e-12,
                "{iso} at {temperature}: {v:?}"
            );
        }
        assert!(matches!(
            p.evaluate(2, 2, 250.0),
            Err(Error::MissingPartitionFunction { .. })
        ));
        assert!(matches!(
            p.evaluate(2, 1, 401.0),
            Err(Error::OutsideTable { .. })
        ));
        assert!(matches!(
            p.evaluate(2, 1, 99.0),
            Err(Error::OutsideTable { .. })
        ));
        assert!(p.insert(2, 3, &[1.0, 1.0], &[1.0, 2.0]).is_err()); // not increasing
        assert!(p.insert(2, 3, &[1.0, 2.0], &[1.0, -2.0]).is_err()); // Q must be > 0
        assert!(p.insert(2, 3, &[1.0], &[1.0]).is_err()); // at least two points
        assert!(p.insert(2, 3, &[1.0, 2.0], &[1.0]).is_err()); // lengths
        assert!(p.insert(2, 3, &[0.0, 2.0], &[1.0, 2.0]).is_err()); // T must be > 0
        let close = f64::from_bits(1e300f64.to_bits() + 1);
        assert!(p.insert(2, 3, &[1e300, close], &[1.0, 2.0]).is_err()); // ln T not increasing
    }
}
