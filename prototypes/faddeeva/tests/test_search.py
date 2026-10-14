import math

import numpy as np
import pytest

from faddeeva_prototype import search
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
