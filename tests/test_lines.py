import json

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
    assert repr(lines.select(isotopologues=1)) == "LineList(1 line, molecules [2], nu 3 to 3 cm^-1)"
    with pytest.raises(TypeError):
        LineList()


def test_molecule_and_isotopologue_numbers_are_checked():
    message = "invalid input: line column `molecule` has values outside 0 to 65535"
    with pytest.raises(absrbrs.InvalidInputError, match=message):
        make(molecule=[2, 2, 70000])
    with pytest.raises(absrbrs.InvalidInputError, match="`isotopologue` has values outside"):
        make(isotopologue=[1, -1, 2])
    with pytest.raises(TypeError):
        make(isotopologue=[1.0, 2.0, 3.0])
    with pytest.raises(
        absrbrs.InvalidInputError, match=r"^invalid input: line column `nu` must be 1-D"
    ):
        make(nu=[[3.0, 1.0, 2.0]])
    with pytest.raises(
        absrbrs.InvalidInputError, match=r"^invalid input: molecule 70000 is outside"
    ):
        make().select(molecule=70000)


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


def test_id_columns_are_read_only():
    lines = make()
    for column in (lines.molecule, lines.isotopologue):
        with pytest.raises(ValueError, match="read-only"):
            column[0] = 1


def test_select_takes_one_isotopologue_and_checks_the_range():
    lines = make()
    np.testing.assert_array_equal(lines.select(isotopologues=2).nu, [1.0])
    assert len(lines.select(wavenumber=np.array([1.5, 3.0]))) == 2
    assert len(lines.select(wavenumber=(np.nan, 3.0))) == 0
    for bad in [("1", "2"), (True, 2.0), (1j, 2.0), ([1.0], 2.0)]:
        with pytest.raises(TypeError, match="wavenumber"):
            lines.select(wavenumber=bad)
    with pytest.raises(ValueError):
        lines.select(wavenumber=(1.0, 2.0, 3.0))
    with pytest.raises(absrbrs.InvalidInputError, match="1-D"):
        lines.select(isotopologues=[[1, 2]])


COLUMNS = dict(
    molecule=np.array([2, 2, 2, 5]),
    isotopologue=np.array([1, 2, 3, 1]),
    nu=np.array([3.0, 1.0, 2.0, 4.0]),
    s_ref=np.array([1e-21, 2e-21, 3e-21, 4e-21]),
    e_lower=np.array([10.0, 20.0, 30.0, 40.0]),
    gamma_air=np.array([0.07, 0.06, 0.05, 0.04]),
    gamma_self=np.array([0.09, 0.08, 0.07, 0.06]),
    n_air=np.array([0.75, 0.7, 0.65, 0.6]),
    delta_air=np.array([-0.004, -0.003, -0.002, -0.001]),
    n_self=np.array([0.6, np.nan, 0.5, 0.4]),
)
FIELDS = ["nu", "s_ref", "e_lower", "gamma_air", "gamma_self", "n_air", "delta_air", "n_self"]

LAYOUTS = {
    "big-endian": lambda a: a.astype(a.dtype.newbyteorder("S")),
    "strided": lambda a: np.repeat(a, 2)[::2],
    "reversed-view": lambda a: a[::-1].copy()[::-1],
}


def assert_same_lines(a, b):
    assert len(a) == len(b)
    np.testing.assert_array_equal(a.molecule, b.molecule)
    np.testing.assert_array_equal(a.isotopologue, b.isotopologue)
    for name in [*FIELDS, "delta_self", "deltap_air"]:
        np.testing.assert_array_equal(getattr(a, name), getattr(b, name))


@pytest.mark.parametrize("layout", LAYOUTS.values(), ids=LAYOUTS.keys())
def test_memory_layout_does_not_matter(layout):
    reference = LineList.from_arrays(**COLUMNS)
    view = {name: layout(values) for name, values in COLUMNS.items()}
    assert not view["nu"].flags.c_contiguous or not view["nu"].dtype.isnative
    assert_same_lines(LineList.from_arrays(**view), reference)
    gamma = np.array([0.01, 0.02, 0.03, 0.04])
    assert_same_lines(
        reference.with_column("gamma_air", layout(gamma)), reference.with_column("gamma_air", gamma)
    )


def test_integer_arguments_are_accepted():
    ints = dict(
        COLUMNS,
        nu=[3, 1, 2, 4],
        e_lower=[10, 20, 30, 40],
        molecule=np.array([2, 2, 2, 5], np.uint8),
    )
    lines = LineList.from_arrays(**ints)
    np.testing.assert_array_equal(lines.nu, [1.0, 2.0, 3.0, 4.0])
    assert lines.nu.dtype == np.float64
    np.testing.assert_array_equal(
        lines.with_column("e_lower", [1, 2, 3, 4]).e_lower, [1.0, 2.0, 3.0, 4.0]
    )


def test_empty_input_gives_empty_outputs():
    empty = LineList.from_arrays(**{name: [] for name in COLUMNS})
    assert len(empty) == 0
    assert empty.nu.shape == empty.n_self.shape == (0,)
    assert empty.nu.dtype == np.float64
    assert empty.molecule.dtype == empty.isotopologue.dtype == np.uint16
    assert len(empty.with_column("gamma_air", [])) == 0
    assert len(empty.select(molecule="CO2", isotopologues=[1], wavenumber=(0.0, 1.0))) == 0


def hapi_header(extra):
    """A HAPI header as hapi.fetch writes it: the fixed part, then extra columns (as in
    absrbrs/tests/io_readers.rs)."""
    fields = [
        ("molec_id", "%2d", 0),
        ("local_iso_id", "%1d", 2),
        ("nu", "%12.6f", 3),
        ("sw", "%10.3E", 15),
        ("a", "%10.3E", 25),
        ("gamma_air", "%5.4f", 35),
        ("gamma_self", "%5.3f", 40),
        ("elower", "%10.4f", 45),
        ("n_air", "%4.2f", 55),
        ("delta_air", "%8.6f", 59),
    ]
    return json.dumps(
        {
            "table_name": "T",
            "order": [name for name, _, _ in fields],
            "format": {name: f for name, f, _ in fields},
            "position": {name: p for name, _, p in fields},
            "extra": extra,
            "extra_format": {name: "%9.6f" for name in extra},
            "extra_separator": ",",
        }
    )


def test_from_hapi(tmp_path):
    (tmp_path / "CO.header").write_text(
        hapi_header(["gamma_sdv_0_air_296", "n_self", "delta_self"])
    )
    a = par_record(5, "1", 2140.0, 1e-19)
    b = par_record(5, "2", 2100.0, 1e-21)
    (tmp_path / "CO.data").write_text(f"{a},0.0500,0.65,-0.001\n{b},#,#,0.003\n")
    lines = LineList.from_hapi(tmp_path, "CO")
    np.testing.assert_array_equal(lines.nu, [2100.0, 2140.0])
    np.testing.assert_array_equal(lines.isotopologue, [2, 1])
    np.testing.assert_array_equal(lines.molecule, [5, 5])
    np.testing.assert_array_equal(lines.n_self, [np.nan, 0.65])  # "#": absent for that line
    np.testing.assert_array_equal(lines.delta_self, [0.003, -0.001])
    assert lines.deltap_air is None
    np.testing.assert_array_equal(lines.gamma_air, [0.07, 0.07])
    np.testing.assert_array_equal(lines.s_ref, [1e-21, 1e-19])
    assert len(LineList.from_hapi(str(tmp_path), "CO")) == 2

    (tmp_path / "X.header").write_text(hapi_header([]))
    with pytest.raises(FileNotFoundError, match=r"X\.data"):
        LineList.from_hapi(tmp_path, "X")
    (tmp_path / "Y.header").write_text("{not json")
    (tmp_path / "Y.data").write_text("")
    with pytest.raises(absrbrs.ParseError):
        LineList.from_hapi(tmp_path, "Y")
