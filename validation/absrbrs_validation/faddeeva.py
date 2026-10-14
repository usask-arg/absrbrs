"""Both Faddeeva tiers against Arb: a dense grid (with every Fast-tier boundary and the double
on its other side) and random points, computed at run time (spec §6.2 metrics)."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from absrbrs import special

TARGETS = {"fast": (1e-6, 1e-5), "accurate": (2e-13, 2e-12)}
X_EDGES = np.array([6.0, 8.0, 9.0, 15.0, 27.0, 50.0, 1750.0])
Y_EDGES = np.array([0.01, 0.1])
RANDOM_POINTS = 200_000
SEED = 20261015


def points() -> NDArray[np.complex128]:
    """The dense grid and the random points, x ≥ 0 and y ≥ 0."""
    xs = np.unique(
        np.concatenate(
            [
                np.arange(385) / 32.0,
                10.0 ** (np.arange(-128, 129) / 16),
                X_EDGES,
                np.nextafter(X_EDGES, -np.inf),
            ]
        )
    )
    ys = np.unique(
        np.concatenate(
            [[0.0], 10.0 ** (np.arange(-180, 37) / 6), Y_EDGES, np.nextafter(Y_EDGES, np.inf)]
        )
    )
    gx, gy = np.meshgrid(xs, ys)
    rng = np.random.default_rng(SEED)
    rx = np.where(
        rng.integers(0, 2, RANDOM_POINTS) == 0,
        rng.uniform(0, 30, RANDOM_POINTS),
        10.0 ** rng.uniform(-8, 8, RANDOM_POINTS),
    )
    ry = np.where(
        rng.integers(0, 10, RANDOM_POINTS) == 0, 0.0, 10.0 ** rng.uniform(-30, 6, RANDOM_POINTS)
    )
    x = np.concatenate([gx.ravel(), rx])
    y = np.concatenate([gy.ravel(), ry])
    z = np.empty(x.shape, dtype=np.complex128)
    z.real, z.imag = x, y
    return z


def measure(z: NDArray[np.complex128], w_ref, dw_ref, tier: str) -> dict:
    """The largest value and derivative errors of ``tier``, where they occur, and how many points
    were exempt (true K below the smallest normal double)."""
    w, dw = special.faddeeva(z, tier=tier)
    k = w_ref.real
    tiny = np.finfo(np.float64).tiny
    with np.errstate(all="ignore"):
        ve = np.where(k < tiny, np.where(np.isfinite(w.real), 0.0, np.inf), np.abs(w.real - k) / k)
        de = np.maximum(np.abs(dw.real - dw_ref.real), np.abs(dw.imag - dw_ref.imag)) / np.abs(
            dw_ref
        )
    ve = np.where(np.isnan(ve), np.inf, ve)
    de = np.where(np.isnan(de), np.inf, de)
    iv, idv = int(np.argmax(ve)), int(np.argmax(de))
    return {
        "tier": tier,
        "points": int(z.size),
        "exempt": int(np.count_nonzero(k < tiny)),
        "max_value_error": float(ve[iv]),
        "at_value": [float(z[iv].real), float(z[iv].imag)],
        "max_derivative_error": float(de[idv]),
        "at_derivative": [float(z[idv].real), float(z[idv].imag)],
        "value_target": TARGETS[tier][0],
        "derivative_target": TARGETS[tier][1],
    }
