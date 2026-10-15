//! Python extension module `absrbrs._absrbrs`.
//!
//! The functions here take and return plain numpy arrays and Python numbers; argument
//! conversion and the user-facing documentation live in the `absrbrs` Python package. Every
//! computing function, and every file read, releases the GIL.

use pyo3::prelude::*;

mod errors;
mod lines;
mod partition;
mod special;
mod testing;

#[pymodule]
fn _absrbrs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add_function(wrap_pyfunction!(special::faddeeva_w, m)?)?;
    m.add_function(wrap_pyfunction!(special::voigt_kernel, m)?)?;
    m.add_class::<lines::PyLineList>()?;
    m.add_function(wrap_pyfunction!(lines::molecule_id, m)?)?;
    m.add_function(wrap_pyfunction!(lines::isotopologue_info, m)?)?;
    m.add_function(wrap_pyfunction!(partition::tips2021, m)?)?;
    testing::register(m)?;
    Ok(())
}
