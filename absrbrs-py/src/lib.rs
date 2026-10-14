//! Python extension module `absrbrs._absrbrs`.
//!
//! The functions here take and return plain numpy arrays and Python numbers; argument
//! conversion and the user-facing documentation live in the `absrbrs` Python package. Every
//! computing function releases the GIL while it computes.

use pyo3::prelude::*;

mod errors;
mod special;
mod testing;

#[pymodule]
fn _absrbrs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add_function(wrap_pyfunction!(special::faddeeva_w, m)?)?;
    m.add_function(wrap_pyfunction!(special::voigt_kernel, m)?)?;
    testing::register(m)?;
    Ok(())
}
