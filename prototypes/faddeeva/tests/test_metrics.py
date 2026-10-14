import numpy as np
import pytest

from faddeeva_prototype.metrics import derivative_error, value_error


def test_value_error_is_relative_to_re_w_and_exempts_subnormal_values():
    w_ref = np.array([2.0 + 5.0j, 1e-310 + 1.0j])
    w = np.array([2.002 + 0.0j, 1.0 + 1.0j])
    err = value_error(w, w_ref)
    assert err[0] == pytest.approx(1e-3)
    assert err[1] == 0.0


def test_derivative_error_is_componentwise_relative_to_the_modulus():
    assert derivative_error(np.array([3.0 + 4.05j]), np.array([3.0 + 4.0j]))[0] == pytest.approx(
        0.01
    )


def test_nan_is_never_hidden():
    assert np.isnan(value_error(np.array([np.nan + 0j]), np.array([1.0 + 0j]))[0])
