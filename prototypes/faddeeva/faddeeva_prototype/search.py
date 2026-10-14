"""Composite (region-switched) algorithms, their errors and costs, and the search over them.

A composite is an ordered list of regions; each point of the (|x|, y) plane uses the method of
the first region that contains it, and the last region contains every point. The errors of
every method at every evaluation point are computed once (`ErrorTable`), so a composite's
errors are a selection from that table.
"""

from __future__ import annotations

import itertools
import math
from dataclasses import asdict, dataclass

import numpy as np
from numpy.typing import NDArray

from faddeeva_prototype.candidates import evaluate
from faddeeva_prototype.metrics import derivative_error, value_error
from faddeeva_prototype.reference import BAND_EDGES, GAUSS_X_MAX

Key = tuple[str, int]

TARGET_VALUE = 1e-6  # spec §6.3, Fast tier
TARGET_DERIVATIVE = 1e-5
#: The prototype must meet half of each target, leaving room for the Rust port's rounding.
MARGIN = 0.5

MAX_DEPTH = 16
#: The method variants the search draws on.
KEYS: tuple[Key, ...] = (
    ("jpole8", 0),
    *(("weideman", n) for n in (24, 32, 40)),
    *(("cf", d) for d in range(MAX_DEPTH + 1)),
    *(("gauss_cf", d) for d in range(MAX_DEPTH + 1)),
    *(("taylor", k) for k in range(1, 6)),
)

#: Cost proxy in flop-equivalents per point, value and derivative together: a complex
#: multiply-add counts 8, a complex division 12, a complex exp 40, Dawson's function 30.
#: It ranks composites; real costs are measured in Rust in phase 1a.
COST_BASE = {"jpole8": 240, "weideman": 56, "cf": 24, "gauss_cf": 94, "taylor": 70}
COST_STEP = {"jpole8": 0, "weideman": 16, "cf": 40, "gauss_cf": 40, "taylor": 24}


def cost(key: Key) -> float:
    """The cost proxy of a method variant."""
    name, param = key
    return COST_BASE[name] + COST_STEP[name] * param


@dataclass(frozen=True)
class Region:
    """Points with y ≤ y_max, x_min ≤ |x| < x_max, r_min ≤ |z| < r_max and |x| + y ≥ s_min."""

    name: str
    method: Key
    y_max: float = math.inf
    x_min: float = 0.0
    x_max: float = math.inf
    r_min: float = 0.0
    r_max: float = math.inf
    s_min: float = 0.0

    def contains(self, z: NDArray[np.complex128]) -> NDArray[np.bool_]:
        x, y, r = np.abs(z.real), z.imag, np.abs(z)
        return (
            (y <= self.y_max)
            & (x >= self.x_min)
            & (x < self.x_max)
            & (r >= self.r_min)
            & (r < self.r_max)
            & (x + y >= self.s_min)
        )


def assign(regions: list[Region], z: NDArray[np.complex128]) -> NDArray[np.intp]:
    """The index of the first region containing each point; every point must be covered."""
    index = np.full(z.shape, -1, dtype=np.intp)
    for i, region in enumerate(regions):
        index[(index < 0) & region.contains(z)] = i
    if np.any(index < 0):
        raise ValueError("the regions do not cover every point; end with a catch-all region")
    return index


class ErrorTable:
    """Value and derivative errors of every method variant at every evaluation point."""

    def __init__(
        self,
        z: NDArray[np.complex128],
        w_ref: NDArray[np.complex128],
        dw_ref: NDArray[np.complex128],
        keys: tuple[Key, ...] = KEYS,
    ) -> None:
        self.z = z
        self.value: dict[Key, NDArray[np.float64]] = {}
        self.derivative: dict[Key, NDArray[np.float64]] = {}
        for key in keys:
            with np.errstate(all="ignore"):
                w, dw = evaluate(key, z)
            self.value[key] = value_error(w, w_ref)
            self.derivative[key] = derivative_error(dw, dw_ref)

    def worst(self, key: Key, mask: NDArray[np.bool_]) -> tuple[float, float]:
        """The largest value and derivative errors of ``key`` over ``mask`` (NaN counts as ∞)."""
        if not mask.any():
            return 0.0, 0.0
        v = self.value[key][mask]
        d = self.derivative[key][mask]
        worst_v = math.inf if np.isnan(v).any() else float(v.max())
        worst_d = math.inf if np.isnan(d).any() else float(d.max())
        return worst_v, worst_d


def meets(errors: tuple[float, float]) -> bool:
    """Whether the errors meet both targets with the margin."""
    return errors[0] <= MARGIN * TARGET_VALUE and errors[1] <= MARGIN * TARGET_DERIVATIVE


def minimal_parameter(
    table: ErrorTable, name: str, mask: NDArray[np.bool_], params: range
) -> int | None:
    """The smallest parameter of method ``name`` that meets the targets on ``mask``; None if
    none does, or if ``mask`` is empty (a region without evaluation points is untested)."""
    if not mask.any():
        return None
    for p in params:
        if meets(table.worst((name, p), mask)):
            return p
    return None


def json_float(v: float) -> float | None:
    """``v``, or None (JSON null) if it is infinite: JSON has no infinity."""
    return v if math.isfinite(v) else None


@dataclass
class Result:
    """A composite with its measured errors and expected cost."""

    regions: list[Region]
    region_errors: list[tuple[float, float]]
    max_value_error: float
    max_derivative_error: float
    expected_cost: float
    occupancy: list[float]

    def to_json(self) -> dict:
        """The composite as strict JSON: unbounded region limits, and errors that are infinite
        (a NaN somewhere), are written as null."""
        return {
            "regions": [
                {
                    **{
                        k: json_float(v) if isinstance(v, float) else v
                        for k, v in asdict(r).items()
                    },
                    "method": list(r.method),
                    "max_value_error": json_float(e[0]),
                    "max_derivative_error": json_float(e[1]),
                    "occupancy": o,
                }
                for r, e, o in zip(self.regions, self.region_errors, self.occupancy, strict=True)
            ],
            "max_value_error": json_float(self.max_value_error),
            "max_derivative_error": json_float(self.max_derivative_error),
            "expected_cost": self.expected_cost,
        }


def holdout_errors(
    regions: list[Region],
    z: NDArray[np.complex128],
    w_ref: NDArray[np.complex128],
    dw_ref: NDArray[np.complex128],
) -> tuple[float, float]:
    """The largest value and derivative errors of a composite at points it was not chosen on
    (NaN counts as ∞)."""
    index = assign(regions, z)
    worst = [0.0, 0.0]
    for i, region in enumerate(regions):
        mask = index == i
        if not mask.any():
            continue
        with np.errstate(all="ignore"):
            w, dw = evaluate(region.method, z[mask])
        for j, err in enumerate((value_error(w, w_ref[mask]), derivative_error(dw, dw_ref[mask]))):
            worst[j] = max(worst[j], math.inf if np.isnan(err).any() else float(err.max()))
    return worst[0], worst[1]


def measure(regions: list[Region], table: ErrorTable, workload: NDArray[np.complex128]) -> Result:
    """Errors over the evaluation points and cost over the workload points of a composite."""
    index = assign(regions, table.z)
    errors = [table.worst(r.method, index == i) for i, r in enumerate(regions)]
    load = assign(regions, workload)
    occupancy = [float(np.mean(load == i)) for i in range(len(regions))]
    return Result(
        regions=regions,
        region_errors=errors,
        max_value_error=max(e[0] for e in errors),
        max_derivative_error=max(e[1] for e in errors),
        expected_cost=sum(o * cost(r.method) for r, o in zip(regions, occupancy, strict=True)),
        occupancy=occupancy,
    )


BASELINES: dict[str, list[Region]] = {
    "jpole8 (sasktran2)": [Region("all", ("jpole8", 0))],
    "hum1wei24 (HAPI)": [
        Region("Humlicek I", ("cf", 1), s_min=15.0),
        Region("Weideman 24", ("weideman", 24)),
    ],
    "Weideman 32": [Region("all", ("weideman", 32))],
    "Weideman 40": [Region("all", ("weideman", 40))],
}


def build(
    table: ErrorTable,
    y_gauss: float,
    x_gauss: float,
    r_core: float,
    core: Key,
    taylor: tuple[float, int] | None,
) -> list[Region] | None:
    """The composite of one template point, with minimal depths; None if a region fails.

    Regions, first match wins: (1) exp(−z²) + continued fraction for y ≤ y_gauss and
    x_gauss ≤ |x| < 27; (2) continued-fraction bands for |z| ≥ r_core; (3) optionally a
    small-y Taylor series of order k for y ≤ y_taylor and |x| < x_gauss; (4) the core method.
    """
    z = table.z
    taken = np.zeros(z.shape, dtype=bool)
    regions: list[Region] = []

    def add(region_without_method: Region, name: str, params: range) -> bool:
        mask = ~taken & region_without_method.contains(z)
        p = minimal_parameter(table, name, mask, params)
        if p is None:
            return False
        regions.append(Region(**{**asdict(region_without_method), "method": (name, p)}))
        taken[mask] = True
        return True

    depths = range(MAX_DEPTH + 1)
    gauss = Region("exp + CF", ("gauss_cf", 0), y_max=y_gauss, x_min=x_gauss, x_max=GAUSS_X_MAX)
    if not add(gauss, "gauss_cf", depths):
        return None
    edges = [r_core, *(e for e in BAND_EDGES if e > r_core), math.inf]
    for lo, hi in itertools.pairwise(edges):
        band = Region(f"CF |z| in [{lo:g}, {hi:g})", ("cf", 0), r_min=lo, r_max=hi)
        if not add(band, "cf", depths):
            return None
    if taylor is not None:
        y_taylor, order = taylor
        small = Region("Taylor", ("taylor", order), y_max=y_taylor, x_max=x_gauss)
        if not add(small, "taylor", range(order, order + 1)):
            return None
    core_region = Region(f"{core[0]} {core[1]}", core)
    if not meets(table.worst(core, ~taken)):
        return None
    regions.append(core_region)
    return regions


GRID = {
    "y_gauss": (0.01, 0.1, 0.5),
    "x_gauss": (4.0, 5.0, 6.0),
    "r_core": (8.0, 10.0, 12.0, 15.0),
    "core": (("weideman", 24), ("weideman", 32), ("weideman", 40)),
    "taylor": (None, (1e-3, 2), (1e-3, 3), (1e-2, 3), (1e-2, 4)),
}


def search(table: ErrorTable, workload: NDArray[np.complex128]) -> list[Result]:
    """Every composite of the template grid that meets the targets, cheapest first."""
    results = []
    for values in itertools.product(*GRID.values()):
        regions = build(table, *values)
        if regions is not None:
            results.append(measure(regions, table, workload))
    return sorted(results, key=lambda r: r.expected_cost)
