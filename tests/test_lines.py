import numpy as np
import pytest

import absrbrs
from absrbrs import LineList


def make(n=3, **extra):
    kw = dict(
        molecule=[2] * n,
        isotopologue=list(range(1, n + 1)),
        nu=[3.0, 1.0, 2.0][:n],
        s_ref=[1e-21, 2e-21, 3e-21][:n],
        e_lower=[10.0] * n,
        gamma_air=[0.07] * n,
        gamma_self=[0.09] * n,
        n_air=[0.75] * n,
        delta_air=[-0.004] * n,
    )
    kw.update(extra)
    return LineList.from_arrays(**kw)


def test_from_arrays_sorts_and_exposes_read_only_columns():
    lines = make()
    assert len(lines) == 3
    np.testing.assert_array_equal(lines.nu, [1.0, 2.0, 3.0])
    np.testing.assert_array_equal(lines.isotopologue, [2, 3, 1])
    assert lines.isotopologue.dtype == np.uint16
    assert lines.n_self is None
    with pytest.raises(ValueError, match="read-only"):
        lines.nu[0] = 5.0


def test_optional_columns_and_with_column():
    lines = make(n_self=[0.6, np.nan, 0.5])
    np.testing.assert_array_equal(lines.n_self, [np.nan, 0.5, 0.6])  # sorted with the rows
    g = lines.with_column("gamma_air", [0.01, 0.02, 0.03])
    np.testing.assert_array_equal(g.gamma_air, [0.01, 0.02, 0.03])
    np.testing.assert_array_equal(lines.gamma_air, [0.07] * 3)  # the original is unchanged
    with pytest.raises(absrbrs.InvalidInputError):
        lines.with_column("gamma_air", [0.01])
    with pytest.raises(absrbrs.InvalidInputError, match="unknown"):
        lines.with_column("sw", [1, 2, 3])


def test_select_by_molecule_name_isotopologues_and_range():
    lines = make()
    assert len(lines.select(molecule="CO2", isotopologues=[1, 3], wavenumber=(1.5, 3.0))) == 2
    assert len(lines.select(molecule=1)) == 0
    with pytest.raises(absrbrs.InvalidInputError, match="molecule"):
        lines.select(molecule="unobtainium")


def test_invalid_arrays():
    with pytest.raises(absrbrs.InvalidInputError):
        make(gamma_air=[0.07, np.nan, 0.07])
    with pytest.raises(TypeError):
        make(nu=[1 + 1j, 2.0, 3.0])


def fixed(v, width, prec):
    """HITRAN's Fortran fixed format: drop the leading zero if needed (".0700", "-.004000")."""
    s = f"{v:.{prec}f}"
    if len(s) > width:
        s = s.replace("0.", ".", 1)
    return s.rjust(width)


def par_record(mol, iso, nu, s):
    r = (
        f"{mol:2d}{iso}{fixed(nu, 12, 6)}{s:10.3E}{1.0:10.3E}{fixed(0.07, 5, 4)}{fixed(0.09, 5, 3)}"
        f"{fixed(100.0, 10, 4)}{fixed(0.75, 4, 2)}{fixed(-0.004, 8, 6)}"
    )
    assert len(r) == 67
    return r.ljust(160)


def test_from_hitran(tmp_path):
    path = tmp_path / "CO2.par"
    path.write_text(
        par_record(2, "1", 6300.5, 1e-23) + "\n" + par_record(2, "A", 6200.0, 1e-25) + "\n"
    )
    lines = LineList.from_hitran(path)
    np.testing.assert_array_equal(lines.isotopologue, [11, 1])
    with pytest.raises(FileNotFoundError):
        LineList.from_hitran(tmp_path / "missing.par")
    (tmp_path / "bad.par").write_text("short\n")
    with pytest.raises(absrbrs.ParseError) as info:
        LineList.from_hitran(tmp_path / "bad.par")
    assert info.value.line_number == 1


def test_select_edge_cases_and_repr():
    lines = make()
    assert len(lines.select(isotopologues=[])) == 0
    assert len(lines.select(wavenumber=(3.0, 1.0))) == 0
    assert repr(lines) == "LineList(3 lines, molecules [2], nu 1 to 3 cm^-1)"
    assert repr(lines.select(molecule=1)) == "LineList(0 lines)"
    with pytest.raises(TypeError):
        LineList()


def test_molecule_and_isotopologue_numbers_are_checked():
    with pytest.raises(absrbrs.InvalidInputError, match="65535"):
        make(molecule=[2, 2, 70000])
    with pytest.raises(absrbrs.InvalidInputError, match="65535"):
        make(isotopologue=[1, -1, 2])
    with pytest.raises(TypeError):
        make(isotopologue=[1.0, 2.0, 3.0])
    with pytest.raises(absrbrs.InvalidInputError, match="1-D"):
        make(nu=[[3.0, 1.0, 2.0]])


def test_isotopologue_metadata():
    assert absrbrs._absrbrs.molecule_id("CO2") == 2
    with pytest.raises(absrbrs.InvalidInputError, match="molecule"):
        absrbrs._absrbrs.molecule_id("co2")
    info = absrbrs._absrbrs.isotopologue_info(2, 1)
    assert info["molecule_name"] == "CO2"
    assert info["global_id"] == 7
    assert info["mass"] == pytest.approx(43.98983)
    assert set(info) == {"formula", "molecule_name", "abundance", "mass", "global_id"}
    with pytest.raises(absrbrs.InvalidInputError, match="isotopologue"):
        absrbrs._absrbrs.isotopologue_info(2, 99)
