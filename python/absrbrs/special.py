"""Special functions: the Faddeeva function and the Voigt kernel (spec §6).

Both accuracy tiers are available: ``"fast"`` (the default; value error ≤ 1e-6 relative to
Re w, derivative error ≤ 1e-5 relative to ``|w′|``) and ``"accurate"`` (≤ 2e-13 and ≤ 2e-12).
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from absrbrs import _absrbrs

__all__ = ["faddeeva", "voigt_kernel"]


def _numeric(name: str, value: ArrayLike) -> NDArray:
    array = np.asarray(value)
    if array.dtype.kind not in "biufc":
        raise TypeError(f"{name} must be numeric, not {array.dtype}")
    return array


def faddeeva(
    z: ArrayLike, *, tier: str = "fast"
) -> tuple[NDArray[np.complex128], NDArray[np.complex128]]:
    """The Faddeeva function w(z) = exp(−z²) erfc(−iz) and its derivative w′(z).

    Parameters
    ----------
    z : array_like
        Points with Im z ≥ 0 (complex or real).
    tier : {"fast", "accurate"}, optional
        The accuracy tier.

    Returns
    -------
    w, dw : numpy.ndarray
        w(z) and w′(z), complex, shaped like ``z``.

    Raises
    ------
    InvalidInputError
        If any Im z < 0, any z is not finite, or ``tier`` is unknown.
    """
    # ravel() copies to a contiguous 1-D array when needed (np.ascontiguousarray would turn a
    # 0-d input into 1-D).
    array = np.asarray(_numeric("z", z), dtype=np.complex128)
    w, dw = _absrbrs.faddeeva_w(array.ravel(), tier)
    return w.reshape(array.shape), dw.reshape(array.shape)


def voigt_kernel(
    x: ArrayLike, y: ArrayLike, *, tier: str = "fast"
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    """The Voigt kernel triplet K = Re w(x + iy), K_x = ∂K/∂x and K_y = ∂K/∂y.

    Parameters
    ----------
    x, y : array_like
        Real points, broadcast together; y ≥ 0.
    tier : {"fast", "accurate"}, optional
        The accuracy tier.

    Returns
    -------
    k, k_x, k_y : numpy.ndarray
        Shaped like the broadcast of ``x`` and ``y``.

    Raises
    ------
    InvalidInputError
        If any y < 0, any value is not finite, or ``tier`` is unknown.
    """
    xb, yb = np.broadcast_arrays(
        _numeric("x", x).astype(np.float64), _numeric("y", y).astype(np.float64)
    )
    shape = xb.shape
    k, k_x, k_y = _absrbrs.voigt_kernel(xb.ravel(), yb.ravel(), tier)
    return k.reshape(shape), k_x.reshape(shape), k_y.reshape(shape)
