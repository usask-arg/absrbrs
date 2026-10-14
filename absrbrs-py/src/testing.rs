//! Hooks for the Python tests (`tests/*.py`), in the submodule `absrbrs._absrbrs._testing`; not
//! part of the absrbrs API.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use crate::errors::to_pyerr;

/// Name under which [`register`] puts the submodule in `sys.modules`.
const NAME: &str = "absrbrs._absrbrs._testing";

/// Adds the `_testing` submodule to `parent` (the extension module).
///
/// It is set as an attribute rather than with `add_submodule`, which would also list it in the
/// parent's `__all__`; and it is put in `sys.modules`, so that `import absrbrs._absrbrs._testing`
/// works too (an extension's submodules are otherwise reachable only as attributes).
pub fn register(parent: &Bound<'_, PyModule>) -> PyResult<()> {
    let py = parent.py();
    let module = PyModule::new(py, NAME)?;
    module.add_function(wrap_pyfunction!(raise_core_error, &module)?)?;
    parent.setattr("_testing", &module)?;
    py.import("sys")?
        .getattr("modules")?
        .set_item(NAME, &module)?;
    Ok(())
}

/// For tests only: raises the exception that the core error `variant` maps to, to test
/// [`to_pyerr`]. Not part of the absrbrs API.
///
/// `variant` names an `absrbrs::Error` variant; the other arguments are its fields (each
/// variant uses the ones it has). `kind` is the `std::io::ErrorKind` of `Io`: `"NotFound"`,
/// `"PermissionDenied"` or `"Other"`.
#[pyfunction]
#[pyo3(signature = (
    variant, message = "", molecule = 0, isotopologue = 0, temperature = 0.0, line_number = 0,
    path = "", kind = "Other"
))]
#[allow(clippy::too_many_arguments)]
pub fn raise_core_error(
    py: Python<'_>,
    variant: &str,
    message: &str,
    molecule: u16,
    isotopologue: u16,
    temperature: f64,
    line_number: usize,
    path: &str,
    kind: &str,
) -> PyResult<()> {
    let err = match variant {
        "InvalidInput" => absrbrs::Error::InvalidInput(message.into()),
        "MissingPartitionFunction" => absrbrs::Error::MissingPartitionFunction {
            molecule,
            isotopologue,
        },
        "OutsideTable" => absrbrs::Error::OutsideTable {
            molecule,
            isotopologue,
            temperature,
        },
        "Parse" => absrbrs::Error::Parse {
            path: path.into(),
            line_number,
            message: message.into(),
        },
        "Io" => absrbrs::Error::Io {
            path: path.into(),
            kind: match kind {
                "NotFound" => std::io::ErrorKind::NotFound,
                "PermissionDenied" => std::io::ErrorKind::PermissionDenied,
                "Other" => std::io::ErrorKind::Other,
                _ => {
                    return Err(PyValueError::new_err(format!(
                        "unknown std::io::ErrorKind {kind:?}"
                    )));
                }
            },
            message: message.into(),
        },
        _ => {
            return Err(PyValueError::new_err(format!(
                "unknown absrbrs::Error variant {variant:?}"
            )));
        }
    };
    Err(to_pyerr(py, err))
}
