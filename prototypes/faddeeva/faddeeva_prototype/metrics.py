"""The error metrics of spec §6.2."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

#: Values whose true K = Re w is below the smallest normal double are exempt (spec §6.2).
TINY = float(np.finfo(np.float64).tiny)


def value_error(w: NDArray[np.complex128], w_ref: NDArray[np.complex128]) -> NDArray[np.float64]:
    """|K − K_ref| / K_ref, with K = Re w; 0 where K_ref < TINY."""
    k = w_ref.real
    with np.errstate(all="ignore"):
        err = np.abs(w.real - k) / k
    return np.where(k >= TINY, err, 0.0)


def derivative_error(
    dw: NDArray[np.complex128], dw_ref: NDArray[np.complex128]
) -> NDArray[np.float64]:
    """max(|ΔRe w′|, |ΔIm w′|) / |w′_ref|: K_x = Re w′ and K_y = −Im w′, each relative to |w′|."""
    with np.errstate(all="ignore"):
        return np.maximum(np.abs(dw.real - dw_ref.real), np.abs(dw.imag - dw_ref.imag)) / np.abs(
            dw_ref
        )
