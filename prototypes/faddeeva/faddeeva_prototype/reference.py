"""Reference values of w and w′: the committed Arb grid, and a denser Arb map cached locally.

The dense map is not committed (it is about 140 000 points); ``dense()`` computes it once with
Arb, in parallel, and caches it in ``prototypes/faddeeva/.cache/dense.npz``.
"""

from __future__ import annotations

import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from flint import acb, arb, ctx
from numpy.typing import NDArray

Complex = NDArray[np.complex128]

PROTOTYPE = Path(__file__).resolve().parents[1]
REPO = PROTOTYPE.parents[1]
COMMITTED = REPO / "reference" / "data" / "faddeeva.jsonl"
CACHE = PROTOTYPE / ".cache" / "dense.npz"

#: Each component used by the metrics is known to this relative accuracy before rounding.
REL_BITS = 60


def complex_points(x: NDArray[np.float64], y: NDArray[np.float64]) -> Complex:
    """x + iy, formed exactly (no arithmetic on the parts)."""
    z = np.empty(np.broadcast(x, y).shape, dtype=np.complex128)
    z.real = x
    z.imag = y
    return z


def load_committed(path: Path = COMMITTED) -> tuple[Complex, Complex, Complex]:
    """The committed reference grid as flat arrays ``(z, w, dw)`` (the nearest doubles)."""
    zs, ws, dws = [], [], []
    for line in path.read_text("utf-8").splitlines():
        if not line.strip():
            continue
        r = json.loads(line)
        x = np.array(r["inputs"]["x"], dtype=np.float64)
        zs.append(complex_points(x, np.full_like(x, r["inputs"]["y"])))
        o = r["outputs"]
        ws.append(complex_points(np.array(o["w_re"]["hi"]), np.array(o["w_im"]["hi"])))
        dws.append(complex_points(np.array(o["dw_re"]["hi"]), np.array(o["dw_im"]["hi"])))
    return np.concatenate(zs), np.concatenate(ws), np.concatenate(dws)


def arb_w(x: float, y: float) -> tuple[float, float, float, float]:
    """(Re w, Im w, Re w′, Im w′) at x + iy, Re w to 2^-60 relative and each part of w′ to
    2^-60 |w′|, computed in Arb at increasing precision."""
    eps = arb(2) ** -REL_BITS
    prec = 128
    while prec <= 16384:
        with ctx.workprec(prec):
            z = acb(x, y)
            w = (-(z * z)).exp() * (-acb(0, 1) * z).erfc()
            dw = -2 * z * w + 2 * acb(0, 1) / arb.pi().sqrt()
            scale = abs(dw)
            ok = (
                w.real.rad() <= abs(w.real.mid()) * eps
                and dw.real.rad() <= scale.mid() * eps
                and dw.imag.rad() <= scale.mid() * eps
            )
            if ok:
                return (
                    float(w.real.mid()),
                    float(w.imag.mid()),
                    float(dw.real.mid()),
                    float(dw.imag.mid()),
                )
        prec *= 2
    raise ArithmeticError(f"w({x!r} + {y!r}i) did not reach 2^-{REL_BITS} by 16384 bits")


def _arb_pair(xy: tuple[float, float]) -> tuple[float, float, float, float]:
    return arb_w(*xy)


def dense_grid() -> Complex:
    """x ∈ {0, 1/32, …, 12} or {10^(k/16) : |k| ≤ 128}; y ∈ {0} or {10^(k/6) : −180 ≤ k ≤ 36}."""
    xs = np.unique(np.concatenate([np.arange(385) / 32.0, 10.0 ** (np.arange(-128, 129) / 16)]))
    ys = np.concatenate([[0.0], 10.0 ** (np.arange(-180, 37) / 6)])
    x, y = np.meshgrid(xs, ys)
    return complex_points(x.ravel(), y.ravel())


def compute(z: Complex, workers: int | None = None) -> tuple[Complex, Complex]:
    """Arb values ``(w, dw)`` at the points ``z``, in parallel."""
    pairs = list(zip(z.real.tolist(), z.imag.tolist(), strict=True))
    with ProcessPoolExecutor(max_workers=workers) as pool:
        out = np.array(list(pool.map(_arb_pair, pairs, chunksize=256)))
    return complex_points(out[:, 0], out[:, 1]), complex_points(out[:, 2], out[:, 3])


def dense(cache: Path = CACHE) -> tuple[Complex, Complex, Complex]:
    """The dense map ``(z, w, dw)``, computed on first use and cached."""
    if cache.exists():
        data = np.load(cache)
        return data["z"], data["w"], data["dw"]
    z = dense_grid()
    w, dw = compute(z)
    cache.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(cache, z=z, w=w, dw=dw)
    return z, w, dw
