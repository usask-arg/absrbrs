"""The line-state reference: tight balls, and identities the formulas must satisfy."""

from fractions import Fraction

from flint import arb, ctx

from absrbrs_reference.generators import line_state


def test_records_cover_the_grid_and_are_tight():
    bits, records = line_state.generate()
    assert bits == 256
    assert len(records) == 4 * 4 * 4 * 3
    assert len({r["case"] for r in records}) == len(records)
    for r in records:
        for out in r["outputs"].values():
            (hi,), (err,) = out["hi"], out["err"]
            assert err <= 2.0**-150 * abs(hi)  # zero only where the value is exactly zero


def test_c2_is_the_exact_si_value():
    assert line_state.C2 == Fraction(662607015 * 299792458 * 100, 1380649 * 10**13)
    assert float(line_state.C2) == 1.4387768775039338


def test_reference_temperature_and_zero_pressure():
    line = dict(zip(line_state.FIELDS, line_state.LINES[1], strict=True))
    with ctx.workprec(256):
        q = line_state.line_shape(296.0, 0.0, 0.5, line, 0.0, 44.0)
        # f = 1 at T_ref with Q(T_ref)/Q(T) = 1, so A·√π·σ_D = 10⁻⁴ S_ref.
        a_check = q["amplitude"] * arb.pi().sqrt() * q["doppler_width"] * 10000
        assert abs(a_check - line_state.exact(line["s_ref"])) < arb(2) ** -200 * line["s_ref"]
        assert q["y"] == 0
        assert q["centre"] == line_state.exact(line["nu"])
