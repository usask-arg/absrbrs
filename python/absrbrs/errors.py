"""Exception hierarchy raised by absrbrs.

Every absrbrs-specific failure is an :class:`AbsrbrsError`:

* :class:`InvalidInputError` (also a :class:`ValueError`): an argument has the right type but
  a value outside its valid range, including the documented size limits.
* :class:`MissingPartitionFunctionError` (also a :class:`LookupError`): no partition function
  is available for an isotopologue in the line list.
* :class:`OutsideTableError` (also a :class:`ValueError`): a temperature is outside the range
  of a partition-function table.
* :class:`ParseError` (also a :class:`ValueError`): a record of a line file is malformed.

A file that cannot be read raises the built-in :class:`OSError`: :class:`FileNotFoundError` if
it does not exist, :class:`PermissionError` if it may not be read.
"""

from __future__ import annotations

__all__ = [
    "AbsrbrsError",
    "InvalidInputError",
    "MissingPartitionFunctionError",
    "OutsideTableError",
    "ParseError",
]


class AbsrbrsError(Exception):
    """Base class for every error raised by absrbrs."""


class InvalidInputError(AbsrbrsError, ValueError):
    """An input is outside its valid range."""


class MissingPartitionFunctionError(AbsrbrsError, LookupError):
    """No partition function is available for an isotopologue.

    Parameters
    ----------
    message : str
        Description of the failure.
    molecule : int, optional
        HITRAN molecule number.
    isotopologue : int, optional
        HITRAN isotopologue number.
    """

    def __init__(
        self, message: str, molecule: int | None = None, isotopologue: int | None = None
    ) -> None:
        super().__init__(message)
        self.molecule = molecule
        self.isotopologue = isotopologue


class OutsideTableError(AbsrbrsError, ValueError):
    """A temperature is outside the range of a partition-function table.

    Parameters
    ----------
    message : str
        Description of the failure.
    molecule : int, optional
        HITRAN molecule number.
    isotopologue : int, optional
        HITRAN isotopologue number.
    temperature : float, optional
        The temperature requested, in K.
    """

    def __init__(
        self,
        message: str,
        molecule: int | None = None,
        isotopologue: int | None = None,
        temperature: float | None = None,
    ) -> None:
        super().__init__(message)
        self.molecule = molecule
        self.isotopologue = isotopologue
        self.temperature = temperature


class ParseError(AbsrbrsError, ValueError):
    """A record of a line file is malformed.

    Parameters
    ----------
    message : str
        Description of the failure.
    path : str, optional
        The file, as given, or a description such as ``<string>`` for input from memory.
    line_number : int, optional
        One-based line number of the record in its file.
    """

    def __init__(
        self, message: str, path: str | None = None, line_number: int | None = None
    ) -> None:
        super().__init__(message)
        self.path = path
        self.line_number = line_number
