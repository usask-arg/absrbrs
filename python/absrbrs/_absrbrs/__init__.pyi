# Stubs for the extension module absrbrs._absrbrs (a stub package, for its _testing submodule).
import os
from typing import Final, TypedDict, final

import numpy as np
from numpy.typing import NDArray

from . import _testing as _testing

__all__ = [
    "LineList",
    "__version__",
    "faddeeva_w",
    "isotopologue_info",
    "molecule_id",
    "tips2021",
    "voigt_kernel",
]

__version__: Final[str]

def faddeeva_w(
    z: NDArray[np.complex128], tier: str
) -> tuple[NDArray[np.complex128], NDArray[np.complex128]]: ...
def voigt_kernel(
    x: NDArray[np.float64], y: NDArray[np.float64], tier: str
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]: ...

@final
class LineList:
    @staticmethod
    def from_hitran(path: str | os.PathLike[str]) -> LineList: ...
    @staticmethod
    def from_hapi(directory: str | os.PathLike[str], table: str) -> LineList: ...
    @staticmethod
    def from_columns(
        molecule: NDArray[np.uint16],
        isotopologue: NDArray[np.uint16],
        nu: NDArray[np.float64],
        s_ref: NDArray[np.float64],
        e_lower: NDArray[np.float64],
        gamma_air: NDArray[np.float64],
        gamma_self: NDArray[np.float64],
        n_air: NDArray[np.float64],
        delta_air: NDArray[np.float64],
        n_self: NDArray[np.float64] | None,
        delta_self: NDArray[np.float64] | None,
        deltap_air: NDArray[np.float64] | None,
    ) -> LineList: ...
    def __len__(self) -> int: ...
    def column(self, name: str) -> NDArray[np.float64] | None: ...
    def molecule(self) -> NDArray[np.uint16]: ...
    def isotopologue(self) -> NDArray[np.uint16]: ...
    def select(
        self,
        molecule: int | None,
        isotopologues: NDArray[np.uint16] | None,
        wavenumber: tuple[float, float] | None,
    ) -> LineList: ...
    def with_column(self, name: str, values: NDArray[np.float64]) -> LineList: ...

class _IsotopologueInfo(TypedDict):
    formula: str
    molecule_name: str
    abundance: float
    mass: float
    global_id: int

def molecule_id(name: str) -> int: ...
def isotopologue_info(molecule: int, isotopologue: int) -> _IsotopologueInfo: ...
def tips2021(
    molecule: int, isotopologue: int, temperature: NDArray[np.float64]
) -> tuple[NDArray[np.float64], NDArray[np.float64]]: ...
