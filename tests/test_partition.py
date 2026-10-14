import numpy as np
import pytest

import absrbrs
from absrbrs import partition


def test_tips2021_values_and_derivative():
    q, d = partition.tips2021("H2O", 1, [296.0, 1.0])
    np.testing.assert_allclose(q, [174.58130, 1.0], rtol=1e-7)
    t = np.array([150.25, 296.0, 777.7])
    _, d = partition.tips2021(2, 1, t)
    qp, _ = partition.tips2021(2, 1, t + 1e-4)
    qm, _ = partition.tips2021(2, 1, t - 1e-4)
    np.testing.assert_allclose(d, (np.log(qp) - np.log(qm)) / 2e-4, rtol=1e-6)


def test_shapes():
    q, d = partition.tips2021(2, 1, np.full((2, 3), 296.0))
    assert q.shape == d.shape == (2, 3)
    q, _ = partition.tips2021(2, 1, 296.0)
    assert q.shape == ()


def test_errors():
    with pytest.raises(absrbrs.OutsideTableError) as info:
        partition.tips2021(3, 1, [296.0, 1500.0])
    assert info.value.temperature == 1500.0
    with pytest.raises(absrbrs.MissingPartitionFunctionError):
        partition.tips2021(10, 3, 296.0)
    with pytest.raises(absrbrs.InvalidInputError, match="molecule"):
        partition.tips2021("XYZ", 1, 296.0)


def test_tables_that_start_above_1_k():
    q, _ = partition.tips2021(31, 2, 30.0)
    assert q > 0
    with pytest.raises(absrbrs.OutsideTableError) as info:
        partition.tips2021(31, 2, 29.0)
    assert (info.value.molecule, info.value.isotopologue) == (31, 2)
    with pytest.raises(absrbrs.OutsideTableError):
        partition.tips2021("CH3", 1, 1.0)
    with pytest.raises(TypeError):
        partition.tips2021(2, 1, [296.0 + 0j])


LAYOUTS = {
    "fortran-order": np.asfortranarray,
    "big-endian": lambda a: a.astype(a.dtype.newbyteorder("S")),
    "reversed-view": lambda a: a[::-1, ::-1],
    "strided": lambda a: np.repeat(a, 2, axis=1)[:, ::2],
}


@pytest.mark.parametrize("layout", LAYOUTS.values(), ids=LAYOUTS.keys())
def test_memory_layout_does_not_matter(layout):
    t = np.linspace(150.0, 900.0, 12).reshape(3, 4)
    view = layout(t)
    for got, expected in zip(
        partition.tips2021(2, 1, view),
        partition.tips2021(2, 1, np.ascontiguousarray(view)),
        strict=True,
    ):
        assert got.shape == (3, 4)
        np.testing.assert_array_equal(got, expected)


def test_integer_arguments_are_accepted():
    for got, expected in zip(
        partition.tips2021(np.int64(2), np.uint16(1), [296, 300]),
        partition.tips2021(2, 1, [296.0, 300.0]),
        strict=True,
    ):
        np.testing.assert_array_equal(got, expected)


def test_empty_input_gives_empty_outputs():
    q, d = partition.tips2021(2, 1, [])
    assert q.shape == d.shape == (0,)
    assert q.dtype == d.dtype == np.float64
    q, d = partition.tips2021(2, 1, np.zeros((0, 3)))
    assert q.shape == d.shape == (0, 3)


def test_argument_types_and_ranges():
    for molecule, isotopologue in [(True, 1), (2, True), (2.0, 1), (2, 1.0)]:
        with pytest.raises(TypeError):
            partition.tips2021(molecule, isotopologue, 296.0)
    with pytest.raises(absrbrs.InvalidInputError, match=r"^invalid input: isotopologue 70000"):
        partition.tips2021(2, 70000, 296.0)
    with pytest.raises(absrbrs.OutsideTableError) as info:
        partition.tips2021(2, 1, np.array([[296.0, 1e5], [1e6, 296.0]]))
    assert info.value.temperature == 1e5  # the first in C order
