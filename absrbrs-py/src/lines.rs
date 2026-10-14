//! Bindings of `absrbrs::lines` and `absrbrs::io`, for `absrbrs.LineList`. Columns are 1-D
//! contiguous arrays (the Python layer converts types and shapes); values are validated by the
//! core. Every column returned is a fresh copy.

use std::path::PathBuf;

use absrbrs::lines::{self, LineColumns, LineField};
use numpy::{IntoPyArray, PyArray1};
use pyo3::prelude::*;
use pyo3::types::PyDict;

use crate::errors::{borrow, copied, invalid, to_pyerr};

/// A copy of the values of `array`.
fn to_vec<T: numpy::Element + Copy>(
    array: &Bound<'_, PyArray1<T>>,
    what: &str,
) -> PyResult<Vec<T>> {
    copied(borrow(array)?.as_slice()?, what)
}

/// The field called `name`, or `InvalidInputError`.
fn field(py: Python<'_>, name: &str) -> PyResult<LineField> {
    LineField::from_name(name).ok_or_else(|| invalid(py, format!("unknown line field {name:?}")))
}

/// A line list: the lines sorted by wavenumber, in structure-of-arrays columns.
#[pyclass(frozen, module = "absrbrs._absrbrs", name = "LineList")]
pub struct PyLineList {
    inner: lines::LineList,
}

impl PyLineList {
    fn wrap(py: Python<'_>, result: absrbrs::Result<lines::LineList>) -> PyResult<PyLineList> {
        result
            .map(|inner| PyLineList { inner })
            .map_err(|e| to_pyerr(py, e))
    }
}

#[pymethods]
impl PyLineList {
    /// The lines of a HITRAN `.par` file (160-byte records).
    #[staticmethod]
    fn from_hitran(py: Python<'_>, path: PathBuf) -> PyResult<PyLineList> {
        let result = py.detach(|| absrbrs::io::read_hitran_par(&path));
        Self::wrap(py, result)
    }

    /// The lines of the HAPI table `table` (`<table>.header` and `<table>.data`) in `directory`.
    #[staticmethod]
    fn from_hapi(py: Python<'_>, directory: PathBuf, table: &str) -> PyResult<PyLineList> {
        let result = py.detach(|| absrbrs::io::read_hapi(&directory, table));
        Self::wrap(py, result)
    }

    /// A list from unsorted columns; the lines are sorted by wavenumber.
    #[staticmethod]
    #[pyo3(signature = (
        molecule, isotopologue, nu, s_ref, e_lower, gamma_air, gamma_self, n_air, delta_air,
        n_self, delta_self, deltap_air
    ))]
    #[allow(clippy::too_many_arguments)]
    fn from_columns(
        py: Python<'_>,
        molecule: Bound<'_, PyArray1<u16>>,
        isotopologue: Bound<'_, PyArray1<u16>>,
        nu: Bound<'_, PyArray1<f64>>,
        s_ref: Bound<'_, PyArray1<f64>>,
        e_lower: Bound<'_, PyArray1<f64>>,
        gamma_air: Bound<'_, PyArray1<f64>>,
        gamma_self: Bound<'_, PyArray1<f64>>,
        n_air: Bound<'_, PyArray1<f64>>,
        delta_air: Bound<'_, PyArray1<f64>>,
        n_self: Option<Bound<'_, PyArray1<f64>>>,
        delta_self: Option<Bound<'_, PyArray1<f64>>>,
        deltap_air: Option<Bound<'_, PyArray1<f64>>>,
    ) -> PyResult<PyLineList> {
        let optional = |array: Option<Bound<'_, PyArray1<f64>>>, what: &str| {
            array.map(|a| to_vec(&a, what)).transpose()
        };
        let columns = LineColumns {
            molecule: to_vec(&molecule, "molecule")?,
            isotopologue: to_vec(&isotopologue, "isotopologue")?,
            nu: to_vec(&nu, "nu")?,
            s_ref: to_vec(&s_ref, "s_ref")?,
            e_lower: to_vec(&e_lower, "e_lower")?,
            gamma_air: to_vec(&gamma_air, "gamma_air")?,
            gamma_self: to_vec(&gamma_self, "gamma_self")?,
            n_air: to_vec(&n_air, "n_air")?,
            delta_air: to_vec(&delta_air, "delta_air")?,
            n_self: optional(n_self, "n_self")?,
            delta_self: optional(delta_self, "delta_self")?,
            deltap_air: optional(deltap_air, "deltap_air")?,
        };
        let result = py.detach(|| lines::LineList::from_columns(columns));
        Self::wrap(py, result)
    }

    fn __len__(&self) -> usize {
        self.inner.len()
    }

    /// A copy of the column `name`, or `None` for an optional column that was not supplied (a
    /// `.par` file never supplies one). An optional column whose entries are all NaN stays an
    /// array.
    fn column<'py>(
        &self,
        py: Python<'py>,
        name: &str,
    ) -> PyResult<Option<Bound<'py, PyArray1<f64>>>> {
        let field = field(py, name)?;
        self.inner
            .column(field)
            .map(|c| Ok(copied(c, name)?.into_pyarray(py)))
            .transpose()
    }

    /// A copy of the HITRAN molecule numbers.
    fn molecule<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray1<u16>>> {
        Ok(copied(self.inner.molecule(), "molecule")?.into_pyarray(py))
    }

    /// A copy of the HITRAN local isotopologue numbers.
    fn isotopologue<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray1<u16>>> {
        Ok(copied(self.inner.isotopologue(), "isotopologue")?.into_pyarray(py))
    }

    /// The lines that match every given filter (`None` places no constraint).
    #[pyo3(signature = (molecule, isotopologues, wavenumber))]
    fn select(
        &self,
        py: Python<'_>,
        molecule: Option<u16>,
        isotopologues: Option<Bound<'_, PyArray1<u16>>>,
        wavenumber: Option<(f64, f64)>,
    ) -> PyResult<PyLineList> {
        let isotopologues = isotopologues.as_ref().map(borrow).transpose()?;
        let isotopologues = isotopologues.as_ref().map(|a| a.as_slice()).transpose()?;
        let inner = py.detach(|| self.inner.select(molecule, isotopologues, wavenumber));
        Ok(PyLineList { inner })
    }

    /// A copy of the list with the column `name` replaced by `values`.
    fn with_column(
        &self,
        py: Python<'_>,
        name: &str,
        values: Bound<'_, PyArray1<f64>>,
    ) -> PyResult<PyLineList> {
        let field = field(py, name)?;
        let values = to_vec(&values, name)?;
        let result = py.detach(|| self.inner.with_column(field, values));
        Self::wrap(py, result)
    }
}

/// The HITRAN molecule number of `name` (exact and case-sensitive, e.g. "CO2").
#[pyfunction]
pub fn molecule_id(py: Python<'_>, name: &str) -> PyResult<u16> {
    lines::molecule_id(name).ok_or_else(|| invalid(py, format!("unknown molecule {name:?}")))
}

/// The properties of an isotopologue: `formula`, `molecule_name`, `abundance`, `mass` (g/mol)
/// and `global_id`. `InvalidInputError` if the table does not have it; a number outside 0 to
/// 65535 raises `OverflowError` (from the argument conversion).
#[pyfunction]
pub fn isotopologue_info(
    py: Python<'_>,
    molecule: u16,
    isotopologue: u16,
) -> PyResult<Bound<'_, PyDict>> {
    let iso = lines::isotopologue(molecule, isotopologue).ok_or_else(|| {
        invalid(
            py,
            format!("unknown isotopologue: molecule {molecule}, isotopologue {isotopologue}"),
        )
    })?;
    let info = PyDict::new(py);
    info.set_item("formula", iso.formula)?;
    info.set_item("molecule_name", iso.molecule_name)?;
    info.set_item("abundance", iso.abundance)?;
    info.set_item("mass", iso.mass)?;
    info.set_item("global_id", iso.global_id)?;
    Ok(info)
}
