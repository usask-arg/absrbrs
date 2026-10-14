//! Bindings of `absrbrs::partition`, for `absrbrs.partition`. The temperatures are a 1-D
//! contiguous array (the Python layer broadcasts and flattens); values are validated by the core.

use absrbrs::partition::{self, PartitionFunction};
use numpy::{IntoPyArray, PyArray1, PyArrayMethods};
use pyo3::exceptions::PyBufferError;
use pyo3::prelude::*;

use crate::errors::{reserved, to_pyerr};

/// TIPS-2021 Q(T) and d ln Q/dT (K⁻¹) of an isotopologue at every temperature (K).
#[pyfunction]
#[allow(clippy::type_complexity)]
pub fn tips2021<'py>(
    py: Python<'py>,
    molecule: u16,
    isotopologue: u16,
    temperature: Bound<'py, PyArray1<f64>>,
) -> PyResult<(Bound<'py, PyArray1<f64>>, Bound<'py, PyArray1<f64>>)> {
    let temperature = temperature
        .try_readonly()
        .map_err(|e| PyBufferError::new_err(e.to_string()))?;
    let temperature = temperature.as_slice()?;
    let mut q = reserved(temperature.len(), "Q")?;
    let mut dlnq_dt = reserved(temperature.len(), "dlnQ/dT")?;
    // The first failure stops the evaluation.
    let result = py.detach(|| {
        let tips = partition::tips2021();
        for &t in temperature {
            let value = tips.evaluate(molecule, isotopologue, t)?;
            q.push(value.q);
            dlnq_dt.push(value.dlnq_dt);
        }
        Ok(())
    });
    result.map_err(|e| to_pyerr(py, e))?;
    Ok((q.into_pyarray(py), dlnq_dt.into_pyarray(py)))
}
