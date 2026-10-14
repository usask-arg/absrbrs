# Stubs for absrbrs._absrbrs._testing: hooks for the tests only, not part of the absrbrs API.
from typing import NoReturn

__all__ = ["raise_core_error"]

def raise_core_error(
    variant: str,
    message: str = "",
    molecule: int = 0,
    isotopologue: int = 0,
    temperature: float = 0.0,
    line_number: int = 0,
    path: str = "",
    kind: str = "Other",
) -> NoReturn: ...
