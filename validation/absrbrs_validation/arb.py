"""Arb values of w(z) and w′(z), in parallel (as in prototypes/faddeeva, without its cache)."""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor

import numpy as np
from flint import acb, arb, ctx
from numpy.typing import NDArray

REL_BITS = 60


def arb_w(xy: tuple[float, float]) -> tuple[float, float, float, float]:
    """(Re w, Im w, Re w′, Im w′) at x + iy.

    Re w is known to 2^-60 relative, and each part of w′ to 2^-60 |w′|.
    """
    x, y = xy
    eps = arb(2) ** -REL_BITS
    bits = 128
    while bits <= 16384:
        with ctx.workprec(bits):
            z = acb(x, y)
            w = (-(z * z)).exp() * (-acb(0, 1) * z).erfc()
            dw = -2 * z * w + 2 * acb(0, 1) / arb.pi().sqrt()
            scale = abs(dw).mid()
            if (
                w.real.rad() <= abs(w.real.mid()) * eps
                and dw.real.rad() <= scale * eps
                and dw.imag.rad() <= scale * eps
            ):
                return (
                    float(w.real.mid()),
                    float(w.imag.mid()),
                    float(dw.real.mid()),
                    float(dw.imag.mid()),
                )
        bits *= 2
    raise ArithmeticError(f"w({x!r} + {y!r}i) did not reach 2^-{REL_BITS} by 16384 bits")


def compute(z: NDArray[np.complex128]) -> tuple[NDArray[np.complex128], NDArray[np.complex128]]:
    """Arb values (w, dw) at every point of z."""
    pairs = list(zip(z.real.tolist(), z.imag.tolist(), strict=True))
    with ProcessPoolExecutor() as pool:
        out = np.array(list(pool.map(arb_w, pairs, chunksize=512)))
    return out[:, 0] + 1j * out[:, 1], out[:, 2] + 1j * out[:, 3]
