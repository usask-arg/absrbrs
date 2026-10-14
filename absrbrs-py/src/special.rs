//! Bindings of `absrbrs::faddeeva`, for `absrbrs.special`. Inputs are 1-D contiguous arrays
//! (the Python layer broadcasts and flattens); values are validated here.

use absrbrs::faddeeva::{self, Complex64, Tier};
use numpy::{IntoPyArray, PyArray1, PyArrayMethods};
use pyo3::exceptions::{PyBufferError, PyMemoryError};
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

/// An empty vector with room for `len` values, or `MemoryError` (spec §9.2: fallible allocation).
fn reserved<T>(len: usize, what: &str) -> PyResult<Vec<T>> {
    let mut v = Vec::new();
    v.try_reserve_exact(len).map_err(|e| {
        PyMemoryError::new_err(format!("cannot allocate {len} values for {what}: {e}"))
    })?;
    Ok(v)
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
    let mut w = reserved(z.len(), "w")?;
    let mut dw = reserved(z.len(), "dw")?;
    // Validation and evaluation share one pass, without the GIL; the first invalid point stops it.
    let first_invalid = py.detach(|| {
        for (i, &z) in z.iter().enumerate() {
            if !(z.im >= 0.0 && z.re.is_finite() && z.im.is_finite()) {
                return Some(i);
            }
            let (value, derivative) = faddeeva::w(z, tier);
            w.push(value);
            dw.push(derivative);
        }
        None
    });
    if let Some(i) = first_invalid {
        let message = format!("z[{i}] = {}: Im z must be >= 0 and z finite", z[i]);
        return Err(invalid(py, message));
    }
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
    let mut k = reserved(x.len(), "K")?;
    let mut k_x = reserved(x.len(), "K_x")?;
    let mut k_y = reserved(x.len(), "K_y")?;
    let first_invalid = py.detach(|| {
        for (i, (&x, &y)) in x.iter().zip(y).enumerate() {
            if !(x.is_finite() && y.is_finite() && y >= 0.0) {
                return Some(i);
            }
            let t = faddeeva::kernel(x, y, tier);
            k.push(t.k);
            k_x.push(t.k_x);
            k_y.push(t.k_y);
        }
        None
    });
    if let Some(i) = first_invalid {
        let message = format!(
            "x[{i}] = {}, y[{i}] = {}: y must be >= 0 and x, y finite",
            x[i], y[i]
        );
        return Err(invalid(py, message));
    }
    Ok((
        k.into_pyarray(py),
        k_x.into_pyarray(py),
        k_y.into_pyarray(py),
    ))
}
