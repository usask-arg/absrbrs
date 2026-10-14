//! Conversion of `absrbrs::Error` into the exception classes defined in `absrbrs.errors`.

use pyo3::exceptions::{
    PyException, PyFileNotFoundError, PyOSError, PyPermissionError, PyRuntimeError,
};
use pyo3::prelude::*;

/// The matching `absrbrs.errors` exception for a core error.
///
/// The message is the error's `Display` text, and the structured fields of a variant become
/// attributes of the exception. `Io` maps to the built-in `FileNotFoundError`,
/// `PermissionError` or (for other kinds of failure) `OSError`. Variants added to the
/// non-exhaustive `absrbrs::Error` later map to `AbsrbrsError` until they get a class of their
/// own. If building the exception fails with an `Exception` (say `absrbrs.errors` cannot be
/// imported), the result is a `RuntimeError` with the core message, caused by that exception;
/// any other failure (`KeyboardInterrupt`, ...) is returned unchanged.
pub fn to_pyerr(py: Python<'_>, err: absrbrs::Error) -> PyErr {
    let message = err.to_string();
    if let absrbrs::Error::Io { kind, .. } = err {
        return match kind {
            std::io::ErrorKind::NotFound => PyFileNotFoundError::new_err(message),
            std::io::ErrorKind::PermissionDenied => PyPermissionError::new_err(message),
            _ => PyOSError::new_err(message),
        };
    }
    let build = || -> PyResult<PyErr> {
        let module = py.import("absrbrs.errors")?;
        let instance = match &err {
            absrbrs::Error::InvalidInput(_) => module
                .getattr("InvalidInputError")?
                .call1((message.as_str(),))?,
            absrbrs::Error::MissingPartitionFunction {
                molecule,
                isotopologue,
            } => module.getattr("MissingPartitionFunctionError")?.call1((
                message.as_str(),
                *molecule,
                *isotopologue,
            ))?,
            absrbrs::Error::OutsideTable {
                molecule,
                isotopologue,
                temperature,
            } => module.getattr("OutsideTableError")?.call1((
                message.as_str(),
                *molecule,
                *isotopologue,
                *temperature,
            ))?,
            absrbrs::Error::Parse {
                path, line_number, ..
            } => module.getattr("ParseError")?.call1((
                message.as_str(),
                path.as_str(),
                *line_number,
            ))?,
            _ => module.getattr("AbsrbrsError")?.call1((message.as_str(),))?,
        };
        Ok(PyErr::from_value(instance))
    };
    build().unwrap_or_else(|failure| {
        if !failure.is_instance_of::<PyException>(py) {
            return failure;
        }
        let e = PyRuntimeError::new_err(message);
        e.set_cause(py, Some(failure));
        e
    })
}
