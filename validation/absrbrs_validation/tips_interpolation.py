"""How TIPS-2021 values between the nodes depend on the interpolation, and how far HAPI's are off.

``uv run --locked --directory validation python -m absrbrs_validation.tips_interpolation``
prints the figures quoted in ``docs/theory/partition-functions.md`` and in the ``Tips2021``
rustdoc. Three interpolations of the same tables are compared:

- absrbrs: cubic Hermite in ln Q on the 1 K grid (``absrbrs.partition.tips2021``);
- ``TIPS_2021.py`` (the Zenodo record's own code): linear in Q on the 1 K grid;
- HAPI 1.3.0.0 (``partitionSum(..., version=2021)``, through ``BD_TIPS_2021_PYTHON`` and
  ``AtoB``): 3- and 4-point Lagrange in Q on its 10 K grid (1, 10, 20, ... K).

At integer temperatures absrbrs returns the stored 1 K values (the source, to f32 rounding), so
HAPI at non-node integer temperatures measures HAPI's interpolation error against the source.
At half-integer temperatures absrbrs is compared with linear interpolation in Q of its own
nodes (what ``TIPS_2021.py`` returns there) and with 6-point Lagrange interpolation of ln Q,
a higher-order interpolant whose spread from absrbrs estimates the uncertainty between nodes.
"""

from __future__ import annotations

import contextlib
import io
import warnings

import numpy as np

import absrbrs
from absrbrs import partition

with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
    warnings.simplefilter("ignore", SyntaxWarning)
    import hapi

# The main atmospheric absorbers (molecules 1-7, isotopologue 1).
MAIN = [(1, 1), (2, 1), (3, 1), (4, 1), (5, 1), (6, 1), (7, 1)]
# Slope discontinuities in the source tables (see the Tips2021 rustdoc).
KINKS = {(47, 1): 650.0, (23, 1): 1000.0, (23, 2): 1000.0}
# Weights of 6-point Lagrange interpolation at the midpoint of the middle interval.
LAGRANGE6 = np.array([3.0, -25.0, 150.0, 150.0, -25.0, 3.0]) / 256.0


def inside(key: tuple[int, int], t: float) -> bool:
    try:
        partition.tips2021(*key, t)
    except absrbrs.OutsideTableError:
        return False
    return True


def tables() -> dict[tuple[int, int], tuple[float, float]]:
    """Every absrbrs TIPS-2021 table with its integer range [t_first, t_last], by probing."""
    out = {}
    for m in range(1, 100):
        for i in range(1, 30):
            try:
                partition.tips2021(m, i, 296.0)
            except absrbrs.MissingPartitionFunctionError:
                continue
            lo = next(float(t) for t in range(1, 297) if inside((m, i), float(t)))
            a, b = 296, 10_000  # inside, outside
            while b - a > 1:
                c = (a + b) // 2
                a, b = (c, b) if inside((m, i), float(c)) else (a, c)
            out[(m, i)] = (lo, float(a))
    return out


def hapi_q(key: tuple[int, int], t: np.ndarray) -> np.ndarray:
    with contextlib.redirect_stdout(io.StringIO()):
        return np.array([hapi.partitionSum(*key, float(x), version=2021) for x in t])


def hapi_vs_source(key, lo, hi):
    """Max |HAPI/source − 1| at the non-node integer temperatures in [lo, hi]."""
    t = np.array([x for x in range(int(lo), int(hi) + 1) if x % 10], dtype=np.float64)
    ours, _ = partition.tips2021(*key, t)
    return float(np.max(np.abs(hapi_q(key, t) / ours - 1.0)))


def midpoints(key, lo, hi):
    """At the half-integer T in [lo, hi]: T, |cubic/linear-in-Q − 1|, |ln cubic − Lagrange6|."""
    t = np.arange(lo, hi + 1.0)
    q, _ = partition.tips2021(*key, t)
    mid = t[:-1] + 0.5
    qm, _ = partition.tips2021(*key, mid)
    linear = np.abs(qm / (0.5 * (q[:-1] + q[1:])) - 1.0)
    lag = np.full(len(mid), np.nan)
    ln_q = np.log(q)
    for k in range(2, len(mid) - 3):
        lag[k] = LAGRANGE6 @ ln_q[k - 2 : k + 4]
    return mid, linear, np.abs(np.log(qm) - lag)


def main() -> None:
    tabs = tables()
    in_hapi = [k for k in sorted(tabs) if k in hapi.TIPS_2021_ISOT_HASH]
    print(f"{len(tabs)} absrbrs tables, {len(in_hapi)} also in HAPI's TIPS-2021")

    print("\nHAPI (version=2021) against the 1 K source, non-node integer T:")
    worst = max((hapi_vs_source(k, 200, 300), k) for k in in_hapi)
    print(f"  200-300 K, all: worst {worst[0]:.2e} at {worst[1]}")
    worst = max((hapi_vs_source(k, 150, 350), k) for k in MAIN)
    print(f"  150-350 K, molecules 1-7 isotopologue 1: worst {worst[0]:.2e} at {worst[1]}")
    low = [k for k in in_hapi if tabs[k][0] <= 2.0]
    worst = max((hapi_vs_source(k, 2, 19), k) for k in low)
    print(f"  2-19 K: worst {worst[0]:.2f} at {worst[1]}")

    print("\nabsrbrs (cubic in ln Q) at half-integer T:")
    stats = {}
    for k, (lo, hi) in sorted(tabs.items()):
        mid, linear, lag = midpoints(k, lo, hi)
        kink = np.abs(mid - KINKS.get(k, -1e9)) < 4.0
        band = (mid >= 150) & (mid <= 350)
        warm = (mid >= 50) & ~kink
        stats[k] = {
            "lin_150_350": linear[band].max(),
            "lin_50": linear[warm].max(),
            "lin_kink": linear[kink].max() if kink.any() else 0.0,
            "lag_50": np.nanmax(lag[warm]),
            "lag_kink": np.nanmax(lag[kink]) if kink.any() else 0.0,
            "lin_low": linear[mid < 20].max() if (mid < 20).any() else 0.0,
            "lag_low": np.nanmax(lag[mid < 20]) if np.isfinite(lag[mid < 20]).any() else 0.0,
        }

    def top(name, fmt=".2e"):
        k = max(stats, key=lambda k: stats[k][name])
        return f"{stats[k][name]:{fmt}} at {k}"

    values = [s["lin_150_350"] for s in stats.values()]
    print(f"  against linear in Q, 150-350 K: median {np.median(values):.2e}, worst", end=" ")
    print(top("lin_150_350"))
    print(f"  against linear in Q, from 50 K, away from the kinks: worst {top('lin_50')}")
    print(f"  against linear in Q, within 4 K of the kinks: worst {top('lin_kink')}")
    values = [s["lag_50"] for s in stats.values()]
    print("  against 6-point Lagrange in ln Q, from 50 K, away from the kinks: median", end=" ")
    print(f"{np.median(values):.1e}, worst {top('lag_50')}")
    print(f"  against 6-point Lagrange in ln Q, within 4 K of the kinks: worst {top('lag_kink')}")
    print(f"  below 20 K, against 6-point Lagrange in ln Q: worst {top('lag_low')}")
    print(f"  below 20 K, against linear in Q: worst {top('lin_low')}", end="")
    rest = max(v["lin_low"] for k, v in stats.items() if k != (57, 1))
    print(f"; {rest:.2f} without 57_1")


if __name__ == "__main__":
    main()
