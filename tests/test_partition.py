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
