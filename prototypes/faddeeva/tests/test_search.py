import json
import math
import re

import numpy as np
import pytest

from faddeeva_prototype import reference, search
from faddeeva_prototype.__main__ import table_markdown
from faddeeva_prototype.reference import arb_w, complex_points
from faddeeva_prototype.search import ErrorTable, Region, assign, cost, minimal_parameter


def test_first_matching_region_wins_and_every_point_must_be_covered():
    regions = [Region("far", ("cf", 1), r_min=10.0), Region("rest", ("weideman", 32))]
    z = np.array([20.0 + 0.0j, 1.0 + 1.0j])
    assert assign(regions, z).tolist() == [0, 1]
    with pytest.raises(ValueError, match="catch-all"):
        assign(regions[:1], z)


def test_regions_use_the_absolute_value_of_x():
    region = Region("r", ("cf", 0), x_min=2.0, x_max=3.0, y_max=1.0)
    z = np.array([-2.5 + 0.5j, 2.5 + 1.5j, 3.0 + 0.0j])
    assert region.contains(z).tolist() == [True, False, False]


def test_minimal_depth_and_costs():
    z = complex_points(np.array([50.0, 2000.0]), np.array([1e-3, 10.0]))
    v = np.array([arb_w(float(q.real), float(q.imag)) for q in z])
    w, dw = complex_points(v[:, 0], v[:, 1]), complex_points(v[:, 2], v[:, 3])
    table = ErrorTable(z, w, dw, keys=tuple(("cf", d) for d in range(5)))
    # The Lorentzian (depth 0) is within half the target at |z| = 2000 but not at |z| = 50.
    assert minimal_parameter(table, "cf", np.array([False, True]), range(5)) == 0
    assert minimal_parameter(table, "cf", np.array([True, False]), range(5)) == 1
    assert cost(("cf", 2)) == 104
    assert cost(("weideman", 32)) == 568


def test_measure_reports_occupancy_that_sums_to_one():
    z = complex_points(np.array([50.0, 2000.0]), np.array([1e-3, 10.0]))
    v = np.array([arb_w(float(q.real), float(q.imag)) for q in z])
    w, dw = complex_points(v[:, 0], v[:, 1]), complex_points(v[:, 2], v[:, 3])
    table = ErrorTable(z, w, dw, keys=(("cf", 1), ("cf", 0)))
    regions = [Region("far", ("cf", 0), r_min=1000.0), Region("rest", ("cf", 1))]
    result = search.measure(regions, table, np.array([3000.0 + 0j, 60.0 + 0j, 70.0 + 0j]))
    assert math.isclose(sum(result.occupancy), 1.0)
    assert result.occupancy == pytest.approx([1 / 3, 2 / 3])
    assert result.expected_cost == pytest.approx((cost(("cf", 0)) + 2 * cost(("cf", 1))) / 3)


def test_holdout_errors_cover_every_region_and_count_nan_as_infinite():
    z = complex_points(np.array([50.0, 2000.0]), np.array([1e-3, 10.0]))
    v = np.array([arb_w(float(q.real), float(q.imag)) for q in z])
    w, dw = complex_points(v[:, 0], v[:, 1]), complex_points(v[:, 2], v[:, 3])
    regions = [Region("far", ("cf", 0), r_min=1000.0), Region("rest", ("cf", 1))]
    value, derivative = search.holdout_errors(regions, z, w, dw)
    assert 0 < value < 5e-7 and 0 < derivative < 5e-7
    # The continued fraction is NaN at z = 0, where w = 1 and w′ = 2i/√π.
    nan_case = search.holdout_errors(
        [Region("all", ("cf", 0))],
        np.array([0j]),
        np.array([1 + 0j]),
        np.array([2j / math.sqrt(math.pi)]),
    )
    assert nan_case == (math.inf, math.inf)


def test_each_band_is_measured_on_its_own_points(monkeypatch):
    # 9 + 0i belongs to the exp + CF region. The continued fraction alone gives K = 0 there, so
    # if the [8, 10) band were measured on all its points, rather than on those no earlier
    # region took, no depth would do.
    monkeypatch.setattr(search, "BAND_EDGES", (8.0, 10.0))
    z = complex_points(np.array([9.0, 9.0, 50.0, 1.0]), np.array([0.0, 1.0, 1.0, 1.0]))
    v = np.array([arb_w(float(q.real), float(q.imag)) for q in z])
    w, dw = complex_points(v[:, 0], v[:, 1]), complex_points(v[:, 2], v[:, 3])
    table = ErrorTable(z, w, dw)
    regions = search.build(table, 0.01, 6.0, 8.0, ("weideman", 24), None)
    assert regions is not None
    assert [r.method[0] for r in regions] == ["gauss_cf", "cf", "cf", "weideman"]
    assert assign(regions, z).tolist() == [0, 1, 2, 3]


def test_a_region_without_points_is_rejected():
    z = np.array([50.0 + 1.0j])
    table = ErrorTable(z, z, z, keys=(("cf", 0),))
    assert minimal_parameter(table, "cf", np.array([False]), range(1)) is None


def test_the_dense_grid_samples_every_threshold_of_the_template():
    assert set(search.GRID["x_gauss"]) | {search.GAUSS_X_MAX} | set(search.BAND_EDGES) <= set(
        reference.X_EDGES
    )
    taylor_y = {t[0] for t in search.GRID["taylor"] if t is not None}
    assert set(search.GRID["y_gauss"]) | taylor_y <= set(reference.Y_EDGES)


def test_report_files_are_strict_json_and_well_formed_markdown():
    z = complex_points(np.array([50.0, 2000.0]), np.array([1e-3, 10.0]))
    v = np.array([arb_w(float(q.real), float(q.imag)) for q in z])
    w, dw = complex_points(v[:, 0], v[:, 1]), complex_points(v[:, 2], v[:, 3])
    table = ErrorTable(z, w, dw, keys=(("cf", 1), ("cf", 0)))
    regions = [Region("CF |z| in [1000, inf)", ("cf", 0), r_min=1000.0), Region("rest", ("cf", 1))]
    result = search.measure(regions, table, z)
    json.dumps(result.to_json(), allow_nan=False)  # raises on an infinite bound
    assert result.to_json()["regions"][0]["r_max"] is None
    markdown = table_markdown(
        {"base": result}, result, {"base": (1.0, 1.0), "**Selected composite**": (1.0, 1.0)}
    )
    for block in markdown.split("\n\n"):
        rows = [line for line in block.splitlines() if line.startswith("|")]
        cells = {len(re.split(r"(?<!\\)\|", row)) for row in rows}
        assert len(cells) <= 1, block
