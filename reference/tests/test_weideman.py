"""The Weideman coefficients: certified, and copied exactly into the Rust Fast tier."""

import re
from pathlib import Path

from flint import ctx

from absrbrs_reference.generators import weideman

RUST = Path(__file__).resolve().parents[2] / "absrbrs" / "src" / "faddeeva" / "fast" / "weideman.rs"


def test_values_are_tight_and_decay():
    with ctx.workprec(weideman.PRECISION_BITS):
        big_l, balls = weideman.coefficients()
        assert all(v.rad() <= 1e-60 for v in [big_l, *balls])  # the Arb balls are tight
    _, records = weideman.generate()
    (r,) = records
    a = r["outputs"]["a"]
    assert len(a["hi"]) == 24
    # err also bounds the remainder of the stored hi + lo + lo2 (about 2^-160 relative).
    assert all(e <= 2.0**-150 * abs(h) for e, h in zip(a["err"], a["hi"], strict=True))
    assert abs(a["hi"][23]) < 1e-9 < abs(a["hi"][0])


def test_rust_constants_are_the_nearest_doubles():
    _, (r,) = weideman.generate()
    text = RUST.read_text("utf-8")
    l_rust = float(
        re.search(r"pub\(crate\) const L: f64 = ([^;]+);", text).group(1).replace("_", "")
    )
    a_block = re.search(r"pub\(crate\) const A: \[f64; N\] = \[(.*?)\];", text, re.S).group(1)
    a_rust = [float(v.replace("_", "")) for v in re.findall(r"-?[0-9][0-9_.e+-]*", a_block)]
    assert l_rust == r["outputs"]["L"]["hi"][0]
    assert a_rust == list(reversed(r["outputs"]["a"]["hi"]))  # highest power first
