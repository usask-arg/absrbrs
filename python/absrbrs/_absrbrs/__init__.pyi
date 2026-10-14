# Stubs for the extension module absrbrs._absrbrs (a stub package, for its _testing submodule).
from typing import Final

import numpy as np
from numpy.typing import NDArray

from . import _testing as _testing

__all__ = ["__version__", "faddeeva_w", "voigt_kernel"]

__version__: Final[str]

def faddeeva_w(
    z: NDArray[np.complex128], tier: str
) -> tuple[NDArray[np.complex128], NDArray[np.complex128]]: ...
def voigt_kernel(
    x: NDArray[np.float64], y: NDArray[np.float64], tier: str
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]: ...
