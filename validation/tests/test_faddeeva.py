"""Both tiers meet their spec §6.3 targets on about 340 000 Arb points computed now."""

import pytest

from absrbrs_validation import arb, faddeeva


@pytest.fixture(scope="session")
def reference():
    z = faddeeva.points()
    w, dw = arb.compute(z)
    return z, w, dw


@pytest.mark.parametrize("tier", ["fast", "accurate"])
def test_tier_meets_its_targets(reference, tier):
    r = faddeeva.measure(*reference, tier)
    assert r["max_value_error"] <= r["value_target"], r
    assert r["max_derivative_error"] <= r["derivative_target"], r


def test_points_cover_the_domain_and_the_boundaries():
    z = faddeeva.points()
    assert z.real.min() == 0.0 and z.real.max() == 1e8
    assert z.imag.min() == 0.0 and z.imag.max() == pytest.approx(1e6)
    for edge in faddeeva.X_EDGES:
        assert (z.real == edge).any()
