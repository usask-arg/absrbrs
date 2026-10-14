"""absrbrs's TIPS-2021 against HAPI's embedded TIPS-2021 (version=2021) at HAPI's own nodes.

HAPI 1.3.0.0 embeds the tables at 10 K steps (1, 10, 20, ... K), where its interpolation returns
the node values, so the comparison is of node values, at HAPI's nodes from 10 K up to 1000 K
(every table, HAPI's and absrbrs's, reaches 1000 K). Nodes below the first temperature of
absrbrs's table (H2S 31_2 and 31_3 start at 30 K) are skipped.

The bound. HAPI's tables hold Q rounded to 7 significant figures (checked below), while absrbrs
stores the 1 K source (6-9 significant figures) as f32. A relative tolerance of 2e-7 cannot hold:
rounding to 7 figures alone moves Q by up to 5e-7 relative (at a mantissa of 1.000000), and every
isotopologue has a node beyond 2e-7 (worst 5.5e-7). So each node is checked against what the
rounding allows and no more: half a unit in the 7th significant figure of HAPI's value, plus the
f32 rounding of absrbrs's value (2^-24 relative), the same form as the converter's cross-check
against the Fortran tables (tools/tips2021.py). For a mantissa near 9.99 this is tighter than
2e-7 (1.1e-7 relative).
"""

import contextlib
import io
import warnings

import numpy as np
import pytest

import absrbrs
from absrbrs import partition

with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
    # HAPI prints a banner on import, and its regular expressions raise SyntaxWarnings.
    warnings.simplefilter("ignore", SyntaxWarning)
    import hapi

T_MAX = 1000.0
# absrbrs has TIPS-2021 tables for these, but HAPI 1.3.0.0's TIPS-2021 does not.
NOT_IN_HAPI = {(56, 1), (57, 1)}
# The tables that start above 10 K, and the HAPI nodes they lose.
LATE_START = {(31, 2): [10.0, 20.0], (31, 3): [10.0, 20.0]}


def has_tips(molecule, isotopologue):
    try:
        partition.tips2021(molecule, isotopologue, 296.0)
    except absrbrs.MissingPartitionFunctionError:
        return False
    return True


def isotopologues():
    """HAPI's isotopologues with an absrbrs TIPS-2021 table and a HAPI TIPS-2021 table."""
    ours = [key for key in sorted(hapi.ISO) if has_tips(*key)]
    return [key for key in ours if key in hapi.TIPS_2021_ISOT_HASH]


def inside(molecule, isotopologue, temperature):
    """Whether `temperature` is inside absrbrs's table, found by probing."""
    try:
        partition.tips2021(molecule, isotopologue, temperature)
    except absrbrs.OutsideTableError:
        return False
    return True


def nodes(molecule, isotopologue):
    """HAPI's nodes in [10 K, 1000 K] with their Q, and the nodes outside absrbrs's table."""
    t = np.asarray(hapi.TIPS_2021_ISOT_HASH[(molecule, isotopologue)], dtype=np.float64)
    q = np.asarray(hapi.TIPS_2021_ISOQ_HASH[(molecule, isotopologue)], dtype=np.float64)
    keep = (t >= 10.0) & (t <= T_MAX)
    t, q = t[keep], q[keep]
    ok = np.array([inside(molecule, isotopologue, float(x)) for x in t])
    return t[ok], q[ok], t[~ok].tolist()


@pytest.mark.parametrize(("molecule", "isotopologue"), isotopologues())
def test_tips_matches_hapi_at_its_nodes(molecule, isotopologue):
    t, table, _ = nodes(molecule, isotopologue)
    assert t[0] == (30.0 if (molecule, isotopologue) in LATE_START else 10.0)
    assert t[-1] == T_MAX
    assert np.array_equal(np.diff(t), np.full(len(t) - 1, 10.0))
    # HAPI's values carry at most 7 significant figures: the premise of the bound.
    assert all(float(f"{v:.6e}") == v for v in table)
    ours, _ = partition.tips2021(molecule, isotopologue, t)
    with contextlib.redirect_stdout(io.StringIO()):
        theirs = np.array(
            [hapi.partitionSum(molecule, isotopologue, float(x), version=2021) for x in t]
        )
    np.testing.assert_array_equal(theirs, table)  # HAPI returns the node values at the nodes
    half_unit = 0.5 * 10.0 ** (np.floor(np.log10(theirs)) - 6)
    bound = half_unit * (1 + 1e-9) + np.abs(ours) * 2.0**-24
    worst = np.argmax(np.abs(ours - theirs) / bound)
    assert np.all(np.abs(ours - theirs) <= bound), (
        f"at {t[worst]} K: absrbrs {ours[worst]!r}, HAPI {theirs[worst]!r}"
    )


def test_coverage():
    cases = isotopologues()
    assert len(cases) == 144
    ours = {key for key in hapi.ISO if has_tips(*key)}
    assert ours - set(cases) == NOT_IN_HAPI
    lost = {key: nodes(*key)[2] for key in cases}
    assert {key: v for key, v in lost.items() if v} == LATE_START
