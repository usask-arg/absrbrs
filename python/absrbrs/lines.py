"""Line lists: spectral lines in HITRAN's parameters, sorted by wavenumber.

Units are HITRAN's: wavenumbers and energies in cm^-1, intensities at 296 K in
cm^-1/(molecule cm^-2), half-widths and shifts at 296 K and 1 atm in cm^-1 atm^-1, and the
temperature coefficient of the air shift in cm^-1 atm^-1 K^-1. Molecules and isotopologues
use HITRAN's numbering (molecule 2 is CO2; local isotopologue 1 is the most abundant).
"""

from __future__ import annotations

import os

import numpy as np
from numpy.typing import ArrayLike, NDArray

from absrbrs import _absrbrs
from absrbrs._convert import molecule as to_molecule
from absrbrs._convert import numbers, read_only, real, reals

__all__ = ["LineList"]


def _column(name: str, values: ArrayLike) -> NDArray[np.float64]:
    return reals(name, f"line column `{name}`", values)


class LineList:
    """Spectral lines, sorted by wavenumber, in structure-of-arrays columns.

    A line list is immutable: :meth:`select` and :meth:`with_column` return new lists. Build
    one with :meth:`from_hitran`, :meth:`from_hapi` or :meth:`from_arrays`; the constructor is
    not public.

    Every column property returns a fresh read-only copy, in the sorted order. The optional
    columns (``n_self``, ``delta_self``, ``deltap_air``) are ``None`` when the column was not
    supplied (a ``.par`` file never supplies them); a NaN entry means "absent for this line",
    and a column of all NaN stays an array. An absent n_self defaults to n_air, an absent
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

        The isotopologue code of a record (column 3) maps to HITRAN's local isotopologue
        number: ``"1"`` to ``"9"`` are 1 to 9, ``"0"`` is 10, and ``"A"``, ``"B"``, ... are 11,
        12, .... Records with S_ref = 0 or nu <= 0 are skipped, and so are empty lines. The
        optional columns are ``None``: a ``.par`` record does not carry them.

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
        return cls._wrap(_absrbrs.LineList.from_hitran(path))

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
        return cls._wrap(_absrbrs.LineList.from_hapi(directory, table))

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
            Intensities at 296 K, cm^-1/(molecule cm^-2) (non-negative). As in HITRAN, they
            include the natural abundance of the isotopologue.
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
            ``None`` means the column is not supplied (absent for every line), and a NaN entry
            absent for that line.

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
        BufferError
            If the extension cannot borrow an array (it is being written to elsewhere).
        MemoryError
            If the list cannot be allocated.
        """

        def optional(name: str, values: ArrayLike | None) -> NDArray[np.float64] | None:
            return None if values is None else _column(name, values)

        inner = _absrbrs.LineList.from_columns(
            numbers("molecule", "line column `molecule`", molecule),
            numbers("isotopologue", "line column `isotopologue`", isotopologue),
            _column("nu", nu),
            _column("s_ref", s_ref),
            _column("e_lower", e_lower),
            _column("gamma_air", gamma_air),
            _column("gamma_self", gamma_self),
            _column("n_air", n_air),
            _column("delta_air", delta_air),
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
        lines = "1 line" if n == 1 else f"{n} lines"
        return f"LineList({lines}, molecules [{molecules}], nu {nu[0]:g} to {nu[-1]:g} cm^-1)"

    def _required(self, name: str) -> NDArray[np.float64]:
        column = self._inner.column(name)
        assert column is not None  # the required columns are always present
        return read_only(column)

    def _optional(self, name: str) -> NDArray[np.float64] | None:
        column = self._inner.column(name)
        return None if column is None else read_only(column)

    @property
    def molecule(self) -> NDArray[np.uint16]:
        """HITRAN molecule numbers (uint16)."""
        return read_only(self._inner.molecule())

    @property
    def isotopologue(self) -> NDArray[np.uint16]:
        """HITRAN local isotopologue numbers (uint16)."""
        return read_only(self._inner.isotopologue())

    @property
    def nu(self) -> NDArray[np.float64]:
        """Line centres, cm^-1, in increasing order."""
        return self._required("nu")

    @property
    def s_ref(self) -> NDArray[np.float64]:
        """Intensities at 296 K, cm^-1/(molecule cm^-2), including the natural abundance."""
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
        """Temperature exponents of the self-broadened half-width, or None if not supplied."""
        return self._optional("n_self")

    @property
    def delta_self(self) -> NDArray[np.float64] | None:
        """Self pressure shifts at 296 K and 1 atm, cm^-1 atm^-1, or None if not supplied."""
        return self._optional("delta_self")

    @property
    def deltap_air(self) -> NDArray[np.float64] | None:
        """Temperature coefficients of the air shift, cm^-1 atm^-1 K^-1, or None if not
        supplied."""
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
        isotopologues : int or array_like of int, optional
            HITRAN local isotopologue numbers (one, or a 1-D sequence); an empty sequence
            selects nothing.
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
            If ``molecule`` is not an int or a str, ``isotopologues`` is not integers, or a
            bound of ``wavenumber`` is not a real number.
        ValueError
            If ``wavenumber`` does not have two items.
        InvalidInputError
            If ``molecule`` is an unknown name, a number does not fit 0 to 65535, or
            ``isotopologues`` has more than one dimension.
        BufferError
            If the extension cannot borrow the isotopologue array.
        MemoryError
            If the result cannot be allocated.
        """
        bounds = None
        if wavenumber is not None:
            lo, hi = wavenumber
            bounds = (real("wavenumber", lo), real("wavenumber", hi))
        isotopes = None
        if isotopologues is not None:
            isotopes = numbers("isotopologues", "isotopologues", np.atleast_1d(isotopologues))
        molecule_id = None if molecule is None else to_molecule(molecule)
        return self._wrap(self._inner.select(molecule_id, isotopes, bounds))

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
            If ``name`` is unknown, or ``values`` is not 1-D or has the wrong length or an
            invalid entry (by the rules of :meth:`from_arrays`).
        BufferError
            If the extension cannot borrow the array (it is being written to elsewhere).
        MemoryError
            If the list cannot be allocated.
        """
        return self._wrap(self._inner.with_column(name, _column(name, values)))
