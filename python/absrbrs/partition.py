"""Partition functions: the built-in TIPS-2021 tables.

TIPS-2021 (Gamache et al. 2021, doi:10.1016/j.jqsrt.2021.107713; MIT licence) tabulates the
total internal partition sum Q(T) of 181 HITRAN isotopologues at 1 K steps, from 1 K to each
table's maximum (1000 to 9000 K). Three tables start higher: H2S isotopologues 2 and 3
(molecule 31) at 30 K, and CH3 (molecule 57) at 2 K. absrbrs stores Q in single precision and
interpolates ln Q in T with cubic Hermite polynomials, so Q and d ln Q/dT are consistent with
each other. From 50 K up (away from kinks in three source tables) the interpolation adds a
median of about 3e-8, and at most 4e-7, to the error of ln Q, and d ln Q/dT carries about 1e-5
relative noise from the data; below about 20 K, where ln Q curves sharply, values between the
nodes are uncertain by 1e-2 to 2e-1.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from absrbrs import _absrbrs
from absrbrs._convert import molecule as to_molecule
from absrbrs._convert import number, numeric

__all__ = ["tips2021"]


def tips2021(
    molecule: int | str, isotopologue: int, temperature: ArrayLike
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """The TIPS-2021 partition function Q(T) and its logarithmic derivative d ln Q/dT.

    Parameters
    ----------
    molecule : int or str
        A HITRAN molecule number, or a molecule name such as ``"H2O"`` (case-sensitive).
    isotopologue : int
        The HITRAN local isotopologue number.
    temperature : array_like
        Temperatures, K, within the isotopologue's table.

    Returns
    -------
    q, dlnq_dt : numpy.ndarray
        Q(T) and d ln Q/dT (K^-1), shaped like ``temperature`` (a 0-d array for a scalar).

    Raises
    ------
    TypeError
        If ``molecule`` is not an int or a str, ``isotopologue`` is not an int (a bool is
        rejected too), or ``temperature`` is not real numbers.
    InvalidInputError
        If ``molecule`` is an unknown name, or a number does not fit 0 to 65535.
    MissingPartitionFunctionError
        If TIPS-2021 has no table for the isotopologue.
    OutsideTableError
        If a temperature (the first in C order) is outside the table, or NaN; its
        ``temperature`` attribute is that temperature.
    BufferError
        If the extension cannot borrow the array (it is being written to elsewhere).
    MemoryError
        If the outputs cannot be allocated.
    """
    array = np.asarray(numeric("temperature", temperature, "iuf", "real numbers"), np.float64)
    q, dlnq_dt = _absrbrs.tips2021(
        to_molecule(molecule), number("isotopologue", isotopologue), array.ravel()
    )
    return q.reshape(array.shape), dlnq_dt.reshape(array.shape)
