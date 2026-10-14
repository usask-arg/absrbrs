"""Argument conversion shared by the absrbrs modules (private).

The converters check types and shapes; values are validated by the extension. Errors raised
here carry the same ``invalid input: `` prefix as those raised by the core.
"""

from __future__ import annotations

import operator
from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray

from absrbrs import _absrbrs
from absrbrs.errors import InvalidInputError

U16_MAX = int(np.iinfo(np.uint16).max)


def invalid_input(message: str) -> InvalidInputError:
    """An InvalidInputError worded like the core's (``invalid input: <message>``)."""
    return InvalidInputError(f"invalid input: {message}")


def numeric(name: str, value: ArrayLike, kinds: str, what: str) -> NDArray[Any]:
    """``value`` as an array, which must have one of the dtype ``kinds`` (TypeError otherwise)."""
    array = np.asarray(value)
    if array.dtype.kind not in kinds:
        raise TypeError(f"{name} must be {what}, not an array of dtype {array.dtype}")
    return array


def real(name: str, value: Any) -> float:
    """A real number (int or float, not bool or str), as a float."""
    array = numeric(name, value, "iuf", "a real number")
    if array.ndim != 0:
        raise TypeError(f"{name} must be a real number, not an array of shape {array.shape}")
    return float(array)


def number(name: str, value: int) -> int:
    """A HITRAN molecule or isotopologue number, checked to fit an unsigned 16-bit integer."""
    if isinstance(value, bool | np.bool_):
        raise TypeError(f"{name} must be an int, not a bool")
    n = operator.index(value)
    if not 0 <= n <= U16_MAX:
        raise invalid_input(f"{name} {n} is outside 0 to {U16_MAX}")
    return n


def molecule(value: int | str) -> int:
    """The HITRAN number of a molecule given by number or by name (``"CO2"``)."""
    if isinstance(value, str):
        return _absrbrs.molecule_id(value)
    return number("molecule", value)


def one_dimensional(what: str, array: NDArray[Any]) -> NDArray[Any]:
    """``array`` as a contiguous 1-D array; ``what`` names it in the error."""
    if array.ndim != 1:
        raise invalid_input(f"{what} must be 1-D, not {array.ndim}-D")
    return np.ascontiguousarray(array)


def numbers(name: str, what: str, values: ArrayLike) -> NDArray[np.uint16]:
    """A contiguous 1-D uint16 array of HITRAN numbers, checked to fit."""
    array = np.asarray(values)
    if array.size == 0 and array.dtype.kind == "f":
        array = array.astype(np.uint16)  # an empty list is float64 to numpy
    array = numeric(name, array, "iu", "integers")
    if array.size and (array.min() < 0 or array.max() > U16_MAX):
        raise invalid_input(f"{what} has values outside 0 to {U16_MAX}")
    return one_dimensional(what, array.astype(np.uint16, copy=False))


def reals(name: str, what: str, values: ArrayLike) -> NDArray[np.float64]:
    """A contiguous 1-D float64 array of real numbers."""
    array = np.asarray(numeric(name, values, "iuf", "real numbers"), dtype=np.float64)
    return one_dimensional(what, array)


def read_only(array: NDArray[Any]) -> NDArray[Any]:
    """``array``, made read-only."""
    array.flags.writeable = False
    return array
