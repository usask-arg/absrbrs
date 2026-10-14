"""Special functions: the Faddeeva function and the Voigt kernel.

The Faddeeva function is w(z) = exp(-z^2) erfc(-iz), and the Voigt kernel is K(x, y) = Re w(x + iy)
with its two partial derivatives. Both functions have two accuracy tiers: ``"fast"`` (the default;
value error <= 1e-6 relative to Re w, derivative error <= 1e-5 relative to ``|w'|``) and
``"accurate"`` (<= 2e-13 and <= 2e-12). The bounds are validated for 0 <= y <= 1e6 and
``|x| <= 1e8``; outside that domain the results are computed without an accuracy guarantee.
"""

from __future__ import annotations

from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray

from absrbrs import _absrbrs

__all__ = ["faddeeva", "voigt_kernel"]


def _numeric(name: str, value: ArrayLike, kinds: str, what: str) -> NDArray[Any]:
    array = np.asarray(value)
    if array.dtype.kind not in kinds:
        raise TypeError(f"{name} must be {what}, not an array of dtype {array.dtype}")
    return array


def faddeeva(
    z: ArrayLike, *, tier: str = "fast"
) -> tuple[NDArray[np.complex128], NDArray[np.complex128]]:
    """The Faddeeva function w(z) = exp(-z^2) erfc(-iz) and its derivative w'(z).

    Parameters
    ----------
    z : array_like
        Points with Im z >= 0 and finite (complex, or real numbers, which are taken as complex
        numbers with Im z = 0).
    tier : {"fast", "accurate"}, optional
        The accuracy tier. The value error bound is relative to Re w, and the derivative error
        bound (componentwise) is relative to ``|w'|``. The bounds are validated for
        0 <= Im z <= 1e6 and ``|Re z| <= 1e8``; outside that domain the results are computed
        without an accuracy guarantee.

    Returns
    -------
    w, dw : numpy.ndarray
        w(z) and w'(z), complex, shaped like ``z`` (a 0-d array for a scalar).

    Raises
    ------
    TypeError
        If ``z`` is not numeric (booleans are rejected too).
    InvalidInputError
        If any Im z < 0, any z is not finite, or ``tier`` is unknown.
    BufferError
        If the extension cannot borrow the array (it is being written to elsewhere).
    MemoryError
        If the outputs cannot be allocated.

    Notes
    -----
    In the Fast tier's Taylor region (Im z <= 0.01 and ``|Re z| < 6``) the returned derivative
    is K_x - iK_y of the (non-analytic) formula used there, which differs from w'(z) by no
    more than the tier's derivative bound.
    """
    # ravel() copies to a contiguous 1-D array when needed (np.ascontiguousarray would turn a
    # 0-d input into 1-D).
    array = np.asarray(_numeric("z", z, "iufc", "numbers"), dtype=np.complex128)
    w, dw = _absrbrs.faddeeva_w(array.ravel(), tier)
    return w.reshape(array.shape), dw.reshape(array.shape)


def voigt_kernel(
    x: ArrayLike, y: ArrayLike, *, tier: str = "fast"
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    """The Voigt kernel triplet K = Re w(x + iy), K_x = dK/dx and K_y = dK/dy.

    Parameters
    ----------
    x, y : array_like
        Real points, broadcast together; y >= 0 and both finite.
    tier : {"fast", "accurate"}, optional
        The accuracy tier. The value error bound is relative to K, and the derivative error
        bound (componentwise, on the pair K_x and K_y) is relative to ``|w'|``. The bounds are
        validated for 0 <= y <= 1e6 and ``|x| <= 1e8``; outside that domain the results are
        computed without an accuracy guarantee.

    Returns
    -------
    k, k_x, k_y : numpy.ndarray
        Shaped like the broadcast of ``x`` and ``y``.

    Raises
    ------
    TypeError
        If ``x`` or ``y`` is not a real number array (complex and boolean are rejected).
    ValueError
        If ``x`` and ``y`` do not broadcast together (raised by numpy).
    InvalidInputError
        If any y < 0, any value is not finite, or ``tier`` is unknown.
    BufferError
        If the extension cannot borrow the array (it is being written to elsewhere).
    MemoryError
        If the outputs cannot be allocated.

    Notes
    -----
    K is the unnormalized Voigt kernel; it integrates to sqrt(pi) over x. It is related to
    ``scipy.special.voigt_profile`` by ``K(x, y) = sigma * sqrt(2 pi) * voigt_profile(u, sigma,
    gamma)`` with ``x = u / (sigma * sqrt(2))`` and ``y = gamma / (sigma * sqrt(2))``.
    """
    xb, yb = np.broadcast_arrays(
        np.asarray(_numeric("x", x, "iuf", "real numbers"), dtype=np.float64),
        np.asarray(_numeric("y", y, "iuf", "real numbers"), dtype=np.float64),
    )
    shape = xb.shape
    k, k_x, k_y = _absrbrs.voigt_kernel(xb.ravel(), yb.ravel(), tier)
    return k.reshape(shape), k_x.reshape(shape), k_y.reshape(shape)
