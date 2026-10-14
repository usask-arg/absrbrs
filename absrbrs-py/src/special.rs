//! Bindings of `absrbrs::faddeeva`, for `absrbrs.special`. Inputs are 1-D contiguous arrays
//! (the Python layer broadcasts and flattens); values are validated here.

use absrbrs::faddeeva::{self, Complex64, Tier};
use numpy::{IntoPyArray, PyArray1, PyArrayMethods};
use pyo3::exceptions::PyBufferError;
use pyo3::prelude::*;

use crate::errors::to_pyerr;

fn tier(py: Python<'_>, name: &str) -> PyResult<Tier> {
    match name {
        "fast" => Ok(Tier::Fast),
        "accurate" => Ok(Tier::Accurate),
        _ => Err(to_pyerr(
            py,
            absrbrs::Error::InvalidInput(format!(
                "tier must be 'fast' or 'accurate', not {name:?}"
            )),
        )),
    }
}

fn invalid(py: Python<'_>, message: String) -> PyErr {
    to_pyerr(py, absrbrs::Error::InvalidInput(message))
}

/// w(z) and w′(z) at every z (Im z ≥ 0, finite).
#[pyfunction]
#[allow(clippy::type_complexity)]
pub fn faddeeva_w<'py>(
    py: Python<'py>,
    z: Bound<'py, PyArray1<Complex64>>,
    tier: &str,
) -> PyResult<(
    Bound<'py, PyArray1<Complex64>>,
    Bound<'py, PyArray1<Complex64>>,
)> {
    let tier = self::tier(py, tier)?;
    let z = z
        .try_readonly()
        .map_err(|e| PyBufferError::new_err(e.to_string()))?;
    let z = z.as_slice()?;
    if let Some(bad) = z
        .iter()
        .find(|z| !(z.im >= 0.0 && z.re.is_finite() && z.im.is_finite()))
    {
        return Err(invalid(
            py,
            format!("Im z must be >= 0 and z finite, got {bad}"),
        ));
    }
    let (w, dw): (Vec<_>, Vec<_>) = py.detach(|| z.iter().map(|&z| faddeeva::w(z, tier)).unzip());
    Ok((w.into_pyarray(py), dw.into_pyarray(py)))
}

/// The kernel triplet (K, K_x, K_y) at every (x, y) pair (y ≥ 0, both finite).
#[pyfunction]
#[allow(clippy::type_complexity)]
pub fn voigt_kernel<'py>(
    py: Python<'py>,
    x: Bound<'py, PyArray1<f64>>,
    y: Bound<'py, PyArray1<f64>>,
    tier: &str,
) -> PyResult<(
    Bound<'py, PyArray1<f64>>,
    Bound<'py, PyArray1<f64>>,
    Bound<'py, PyArray1<f64>>,
)> {
    let tier = self::tier(py, tier)?;
    let x = x
        .try_readonly()
        .map_err(|e| PyBufferError::new_err(e.to_string()))?;
    let y = y
        .try_readonly()
        .map_err(|e| PyBufferError::new_err(e.to_string()))?;
    let (x, y) = (x.as_slice()?, y.as_slice()?);
    if x.len() != y.len() {
        return Err(invalid(
            py,
            format!("x and y have {} and {} values", x.len(), y.len()),
        ));
    }
    if let Some((a, b)) = x
        .iter()
        .zip(y)
        .find(|(a, b)| !(a.is_finite() && b.is_finite() && **b >= 0.0))
    {
        return Err(invalid(
            py,
            format!("y must be >= 0 and x, y finite, got x = {a}, y = {b}"),
        ));
    }
    let (k, (k_x, k_y)): (Vec<_>, (Vec<_>, Vec<_>)) = py.detach(|| {
        x.iter()
            .zip(y)
            .map(|(&x, &y)| {
                let t = faddeeva::kernel(x, y, tier);
                (t.k, (t.k_x, t.k_y))
            })
            .unzip()
    });
    Ok((
        k.into_pyarray(py),
        k_x.into_pyarray(py),
        k_y.into_pyarray(py),
    ))
}
