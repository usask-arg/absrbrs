"""Both tiers meet their spec §6.3 targets on about 350 000 Arb points computed now."""

import json

import numpy as np
import pytest

from absrbrs import special
from absrbrs_validation import __main__ as cli
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
    x, y = z.real, z.imag
    for edge in faddeeva.X_EDGES:
        assert (x == edge).any()
        assert (x == np.nextafter(edge, -np.inf)).any()
    for edge in faddeeva.Y_EDGES:
        assert (y == edge).any()
        assert (y == np.nextafter(edge, np.inf)).any()
    for edge in faddeeva.JOHNSON_X:
        assert (x == edge).any()
        assert (x == np.nextafter(edge, -np.inf)).any()
    for edge in faddeeva.JOHNSON_Y:
        assert (y == edge).any()
        assert (y == np.nextafter(edge, np.inf)).any()
    assert (y == 0.0).any()
    assert (np.isclose(y, 1e-30, rtol=1e-6) & (y > 0)).any()


def test_cli_writes_the_report_and_exits_zero(monkeypatch, tmp_path):
    z = np.array([0.5 + 0.5j, 3 + 0.1j, 20 + 1j, 100 + 2j, 1 + 0j], dtype=np.complex128)
    monkeypatch.setattr(faddeeva, "points", lambda: z)
    monkeypatch.setattr(arb, "compute", lambda zz: special.faddeeva(zz, tier="accurate"))
    monkeypatch.setattr(cli, "DOCS", tmp_path)
    assert cli.main(["faddeeva"]) == 0
    assert (tmp_path / "faddeeva-accuracy-table.md").read_text().startswith("| Tier")
    report = json.loads((tmp_path / "data" / "faddeeva-accuracy.json").read_text())
    assert [r["tier"] for r in report] == ["fast", "accurate"]
