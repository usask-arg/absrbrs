"""Line lists: spectral lines in HITRAN's parameters, sorted by wavenumber.

Units are HITRAN's: wavenumbers and energies in cm^-1, intensities at 296 K in
cm^-1/(molecule cm^-2), half-widths and shifts at 296 K and 1 atm in cm^-1 atm^-1, and the
temperature coefficient of the air shift in cm^-1 atm^-1 K^-1. Molecules and isotopologues
use HITRAN's numbering (molecule 2 is CO2; local isotopologue 1 is the most abundant).
"""

from __future__ import annotations

import operator
import os
from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray

from absrbrs import _absrbrs
from absrbrs.errors import InvalidInputError
from absrbrs.special import _numeric

__all__ = ["LineList"]

_U16_MAX = np.iinfo(np.uint16).max


def _number(name: str, value: int) -> int:
    """A HITRAN molecule or isotopologue number, checked to fit an unsigned 16-bit integer."""
    if isinstance(value, bool | np.bool_):
        raise TypeError(f"{name} must be an int, not a bool")
    number = operator.index(value)
    if not 0 <= number <= _U16_MAX:
        raise InvalidInputError(f"{name} {number} is outside 0 to {_U16_MAX}")
    return number


def _molecule(molecule: int | str) -> int:
    """The HITRAN number of a molecule given by number or by name (``"CO2"``)."""
    if isinstance(molecule, str):
        return _absrbrs.molecule_id(molecule)
    return _number("molecule", molecule)


def _one_dimensional(name: str, array: NDArray[Any]) -> NDArray[Any]:
    if array.ndim != 1:
        raise InvalidInputError(f"{name} must be 1-D, not {array.ndim}-D")
    return np.ascontiguousarray(array)


def _numbers(name: str, values: ArrayLike) -> NDArray[np.uint16]:
    """A 1-D uint16 array of HITRAN numbers, checked to fit."""
    array = np.asarray(values)
    if array.size == 0 and array.dtype.kind == "f":
        array = array.astype(np.uint16)  # an empty list is float64 to numpy
    array = _numeric(name, array, "iu", "integers")
    if array.size and (array.min() < 0 or array.max() > _U16_MAX):
        raise InvalidInputError(f"{name} has values outside 0 to {_U16_MAX}")
    return _one_dimensional(name, array.astype(np.uint16))


def _reals(name: str, values: ArrayLike) -> NDArray[np.float64]:
    """A 1-D float64 array of real numbers."""
    array = np.asarray(_numeric(name, values, "iuf", "real numbers"), dtype=np.float64)
    return _one_dimensional(name, array)


def _read_only(array: NDArray[Any]) -> NDArray[Any]:
    array.flags.writeable = False
    return array


class LineList:
    """Spectral lines, sorted by wavenumber, in structure-of-arrays columns.

    A line list is immutable: :meth:`select` and :meth:`with_column` return new lists. Build
    one with :meth:`from_hitran`, :meth:`from_hapi` or :meth:`from_arrays`; the constructor is
    not public.

    Every column property returns a fresh read-only copy, in the sorted order. The optional
    columns (``n_self``, ``delta_self``, ``deltap_air``) are ``None`` when absent for every line,
    and a NaN entry means "absent for this line"; an absent n_self defaults to n_air, an absent
    delta_self to delta_air, and an absent deltap_air to 0.
    """

    __slots__ = ("_inner",)

    _inner: _absrbrs.LineList

    def __init__(self) -> None:
        raise TypeError(
            "use LineList.from_hitran, LineList.from_hapi or LineList.from_arrays to build a "
            "line list"
        )

    @classmethod
    def _wrap(cls, inner: _absrbrs.LineList) -> LineList:
        lines = object.__new__(cls)
        lines._inner = inner
        return lines

    @classmethod
    def from_hitran(cls, path: str | os.PathLike[str]) -> LineList:
        """Read a HITRAN ``.par`` file (160-byte fixed-width records).

        Records with S_ref = 0 or nu <= 0 are skipped, and so are empty lines. The optional
        columns are absent: a ``.par`` record does not carry them.

        Parameters
        ----------
        path : str or os.PathLike
            The file.

        Returns
        -------
        LineList

        Raises
        ------
        FileNotFoundError, PermissionError, OSError
            If the file cannot be read.
        ParseError
            If a record is malformed; ``line_number`` is its line in the file.
        """
        return cls._wrap(_absrbrs.LineList.from_hitran(os.fspath(path)))

    @classmethod
    def from_hapi(cls, directory: str | os.PathLike[str], table: str) -> LineList:
        """Read a HAPI table: ``<table>.header`` and ``<table>.data`` in ``directory``.

        The fixed-width fields are located from the header. n_self, delta_self and deltap_air
        are read when the header lists them, in ``order`` (as ``hapi.cache2storage`` writes
        them) or in ``extra`` (as ``hapi.fetch`` writes them); a missing value (``#``) is NaN,
        absent for that line. Records with S_ref = 0 or nu <= 0 are skipped.

        Parameters
        ----------
        directory : str or os.PathLike
            The directory holding the table's two files.
        table : str
            The table name.

        Returns
        -------
        LineList

        Raises
        ------
        FileNotFoundError, PermissionError, OSError
            If either file cannot be read.
        ParseError
            If the header or a record is malformed.
        """
        return cls._wrap(_absrbrs.LineList.from_hapi(os.fspath(directory), table))

    @classmethod
    def from_arrays(
        cls,
        *,
        molecule: ArrayLike,
        isotopologue: ArrayLike,
        nu: ArrayLike,
        s_ref: ArrayLike,
        e_lower: ArrayLike,
        gamma_air: ArrayLike,
        gamma_self: ArrayLike,
        n_air: ArrayLike,
        delta_air: ArrayLike,
        n_self: ArrayLike | None = None,
        delta_self: ArrayLike | None = None,
        deltap_air: ArrayLike | None = None,
    ) -> LineList:
        """Build a line list from 1-D columns, in any order; the lines are sorted by nu.

        Parameters
        ----------
        molecule, isotopologue : array_like of int
            HITRAN molecule and local isotopologue numbers.
        nu : array_like
            Line centres, cm^-1 (positive).
        s_ref : array_like
            Intensities at 296 K, cm^-1/(molecule cm^-2) (non-negative).
        e_lower : array_like
            Lower-state energies, cm^-1.
        gamma_air, gamma_self : array_like
            Air- and self-broadened half-widths at 296 K and 1 atm, cm^-1 atm^-1
            (non-negative).
        n_air : array_like
            Temperature exponents of the air-broadened half-width.
        delta_air : array_like
            Air pressure shifts at 296 K and 1 atm, cm^-1 atm^-1.
        n_self, delta_self, deltap_air : array_like, optional
            Temperature exponent of the self-broadened half-width, self pressure shift
            (cm^-1 atm^-1) and temperature coefficient of the air shift (cm^-1 atm^-1 K^-1).
            ``None`` means absent for every line, and a NaN entry absent for that line.

        Returns
        -------
        LineList

        Raises
        ------
        TypeError
            If a column is not real numbers (``molecule`` and ``isotopologue``: integers).
        InvalidInputError
            If a column is not 1-D, the columns differ in length, a required value is not
            finite, nu <= 0, S_ref or a half-width is negative, an optional value is infinite,
            or a molecule or isotopologue number does not fit 0 to 65535.
        """

        def optional(name: str, values: ArrayLike | None) -> NDArray[np.float64] | None:
            return None if values is None else _reals(name, values)

        inner = _absrbrs.LineList.from_columns(
            _numbers("molecule", molecule),
            _numbers("isotopologue", isotopologue),
            _reals("nu", nu),
            _reals("s_ref", s_ref),
            _reals("e_lower", e_lower),
            _reals("gamma_air", gamma_air),
            _reals("gamma_self", gamma_self),
            _reals("n_air", n_air),
            _reals("delta_air", delta_air),
            optional("n_self", n_self),
            optional("delta_self", delta_self),
            optional("deltap_air", deltap_air),
        )
        return cls._wrap(inner)

    def __len__(self) -> int:
        return len(self._inner)

    def __repr__(self) -> str:
        n = len(self)
        if n == 0:
            return "LineList(0 lines)"
        molecules = ", ".join(str(m) for m in np.unique(self.molecule))
        nu = self.nu
        return f"LineList({n} lines, molecules [{molecules}], nu {nu[0]:g} to {nu[-1]:g} cm^-1)"

    def _required(self, name: str) -> NDArray[np.float64]:
        column = self._inner.column(name)
        assert column is not None  # the required columns are always present
        return _read_only(column)

    def _optional(self, name: str) -> NDArray[np.float64] | None:
        column = self._inner.column(name)
        return None if column is None else _read_only(column)

    @property
    def molecule(self) -> NDArray[np.uint16]:
        """HITRAN molecule numbers (uint16)."""
        return _read_only(self._inner.molecule())

    @property
    def isotopologue(self) -> NDArray[np.uint16]:
        """HITRAN local isotopologue numbers (uint16)."""
        return _read_only(self._inner.isotopologue())

    @property
    def nu(self) -> NDArray[np.float64]:
        """Line centres, cm^-1, in increasing order."""
        return self._required("nu")

    @property
    def s_ref(self) -> NDArray[np.float64]:
        """Intensities at 296 K, cm^-1/(molecule cm^-2)."""
        return self._required("s_ref")

    @property
    def e_lower(self) -> NDArray[np.float64]:
        """Lower-state energies, cm^-1."""
        return self._required("e_lower")

    @property
    def gamma_air(self) -> NDArray[np.float64]:
        """Air-broadened half-widths at 296 K and 1 atm, cm^-1 atm^-1."""
        return self._required("gamma_air")

    @property
    def gamma_self(self) -> NDArray[np.float64]:
        """Self-broadened half-widths at 296 K and 1 atm, cm^-1 atm^-1."""
        return self._required("gamma_self")

    @property
    def n_air(self) -> NDArray[np.float64]:
        """Temperature exponents of the air-broadened half-width."""
        return self._required("n_air")

    @property
    def delta_air(self) -> NDArray[np.float64]:
        """Air pressure shifts at 296 K and 1 atm, cm^-1 atm^-1."""
        return self._required("delta_air")

    @property
    def n_self(self) -> NDArray[np.float64] | None:
        """Temperature exponents of the self-broadened half-width, or None if absent."""
        return self._optional("n_self")

    @property
    def delta_self(self) -> NDArray[np.float64] | None:
        """Self pressure shifts at 296 K and 1 atm, cm^-1 atm^-1, or None if absent."""
        return self._optional("delta_self")

    @property
    def deltap_air(self) -> NDArray[np.float64] | None:
        """Temperature coefficients of the air shift, cm^-1 atm^-1 K^-1, or None if absent."""
        return self._optional("deltap_air")

    def select(
        self,
        *,
        molecule: int | str | None = None,
        isotopologues: ArrayLike | None = None,
        wavenumber: tuple[float, float] | None = None,
    ) -> LineList:
        """The lines that match every given filter.

        Parameters
        ----------
        molecule : int or str, optional
            A HITRAN molecule number, or a molecule name such as ``"CO2"`` (case-sensitive).
        isotopologues : array_like of int, optional
            HITRAN local isotopologue numbers; an empty sequence selects nothing.
        wavenumber : tuple of float, optional
            An inclusive range ``(lo, hi)`` of line centres, cm^-1; ``lo > hi`` or a NaN bound
            selects nothing.

        Returns
        -------
        LineList
            The matching lines, still sorted. ``None`` filters place no constraint.

        Raises
        ------
        TypeError
            If ``isotopologues`` is not integers.
        InvalidInputError
            If ``molecule`` is an unknown name, or a number does not fit 0 to 65535.
        """
        lo_hi = None if wavenumber is None else (float(wavenumber[0]), float(wavenumber[1]))
        inner = self._inner.select(
            None if molecule is None else _molecule(molecule),
            None
            if isotopologues is None
            else _numbers("isotopologues", np.atleast_1d(np.asarray(isotopologues))),
            lo_hi,
        )
        return self._wrap(inner)

    def with_column(self, name: str, values: ArrayLike) -> LineList:
        """A copy of the list with one column replaced.

        Replacing ``"nu"`` re-sorts the lines, so arrays taken from this list in the old order
        no longer line up with the new one. Replacing an optional column makes it present.

        Parameters
        ----------
        name : str
            The column: ``"nu"``, ``"s_ref"``, ``"e_lower"``, ``"gamma_air"``,
            ``"gamma_self"``, ``"n_air"``, ``"n_self"``, ``"delta_air"``, ``"deltap_air"`` or
            ``"delta_self"``.
        values : array_like
            The new values, one per line in the current (sorted) order.

        Returns
        -------
        LineList

        Raises
        ------
        TypeError
            If ``values`` is not real numbers.
        InvalidInputError
            If ``name`` is unknown, or ``values`` has the wrong length or an invalid entry (by
            the rules of :meth:`from_arrays`).
        """
        return self._wrap(self._inner.with_column(name, _reals(name, values)))
