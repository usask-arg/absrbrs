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
    assert np.isnan(derivative_error(np.array([np.nan + 1j]), np.array([1.0 + 1j]))[0])


def test_non_finite_values_count_even_at_exempt_points():
    w_ref = np.full(3, 1e-320 + 1j)
    err = value_error(np.array([np.nan + 0j, np.inf + 0j, 1.0 + 0j]), w_ref)
    assert np.isnan(err[0]) and np.isnan(err[1]) and err[2] == 0.0


def test_a_nan_reference_is_never_exempt():
    assert np.isnan(value_error(np.array([1.0 + 0j]), np.array([complex(np.nan, 1.0)]))[0])
