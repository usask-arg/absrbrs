"""Each candidate against Arb where it is meant to be accurate, and its derivative against
finite differences of the candidate itself (the derivative must be that of the formula)."""

import numpy as np
import pytest

from faddeeva_prototype import candidates
from faddeeva_prototype.metrics import derivative_error, value_error
from faddeeva_prototype.reference import arb_w, complex_points

INV_SQRT_PI = 1 / np.sqrt(np.pi)


def arb_values(z):
    v = np.array([arb_w(float(q.real), float(q.imag)) for q in np.ravel(z)])
    return complex_points(v[:, 0], v[:, 1]), complex_points(v[:, 2], v[:, 3])


def grid(xs, ys):
    x, y = np.meshgrid(xs, ys)
    return complex_points(x.ravel(), y.ravel())


def relative(a, b):
    return np.max(np.abs(a - b) / np.abs(b))


@pytest.mark.parametrize(("n", "tol_w", "tol_dw"), [(32, 1e-12, 1e-10), (40, 1e-14, 1e-13)])
def test_weideman_matches_arb_in_modulus(n, tol_w, tol_dw):
    # |Δw|/|w| and |Δw′|/|w′| for 0 ≤ x ≤ 8, 1e-6 ≤ y ≤ 10; the spec metric (relative to K) is
    # what the search measures.
    z = grid(np.linspace(0, 8, 33), np.logspace(-6, 1, 15))
    w_ref, dw_ref = arb_values(z)
    w, dw = candidates.weideman(z, n)
    assert relative(w, w_ref) < tol_w
    assert relative(dw, dw_ref) < tol_dw


def test_continued_fraction_converges_for_large_z():
    angles = np.linspace(0, np.pi / 2, 9)
    polar = np.concatenate([r * np.exp(1j * angles) for r in (20.0, 40.0, 100.0)])
    z = complex_points(np.abs(polar.real), np.abs(polar.imag))
    w_ref, dw_ref = arb_values(z)
    w, dw = candidates.laplace_cf(z, 8)
    assert relative(w, w_ref) < 1e-14
    assert relative(dw, dw_ref) < 1e-14


def test_shallow_continued_fractions_are_the_classical_formulas():
    z = np.array([3.0 + 2.0j, 50.0 + 0.1j])
    w0, dw0 = candidates.laplace_cf(z, 0)
    np.testing.assert_allclose(w0, 1j * INV_SQRT_PI / z, rtol=1e-15)
    np.testing.assert_allclose(dw0, -1j * INV_SQRT_PI / z**2, rtol=1e-15)
    w1, _ = candidates.laplace_cf(z, 1)
    np.testing.assert_allclose(w1, 1j * z * INV_SQRT_PI / (z**2 - 0.5), rtol=1e-14)


def test_gaussian_plus_continued_fraction_near_the_real_axis():
    z = grid(np.linspace(7, 12, 11), [0.0, 1e-20, 1e-8, 1e-2])
    w_ref, dw_ref = arb_values(z)
    w, dw = candidates.gauss_cf(z, 8)
    assert value_error(w, w_ref).max() < 1e-9
    assert derivative_error(dw, dw_ref).max() < 1e-9


@pytest.mark.parametrize(("y", "tol_w", "tol_dw"), [(1e-4, 1e-13, 1e-11), (1e-2, 2e-7, 5e-6)])
def test_small_y_taylor_series(y, tol_w, tol_dw):
    # Order 3 over the region the composites may give it (|x| < 6, y ≤ 0.01). The derivative is
    # that of the formula: its y-part is a series one order shorter, hence the looser tolerance.
    z = complex_points(np.linspace(0, 6, 61), np.full(61, y))
    w_ref, dw_ref = arb_values(z)
    w, dw = candidates.small_y_taylor(z, 3)
    assert value_error(w, w_ref).max() < tol_w
    assert derivative_error(dw, dw_ref).max() < tol_dw


def test_jpole8_reproduces_the_sasktran2_error():
    # Spec Appendix A: about 1.3e-3 relative near x = 9.76 at y = 1e-3.
    z = complex_points(np.linspace(0, 20, 2001), np.full(2001, 1e-3))
    w_ref, _ = arb_values(z)
    w, _ = candidates.jpole8(z)
    err = value_error(w, w_ref)
    assert 1e-3 < err.max() < 2e-3
    assert abs(z[err.argmax()].real - 9.76) < 0.05


FORMULA_POINTS = [
    (("weideman", 32), 1.3 + 0.7j),
    (("cf", 4), 25.0 + 3.0j),
    (("gauss_cf", 6), 6.0 + 0.01j),
    (("taylor", 3), 2.0 + 1e-3j),
    (("taylor", 3), 2.6 + 1e-2j),
    (("jpole8", 0), 1.3 + 0.7j),
]


def k_x_and_k_y(key, z0):
    """The formula's K_x and K_y by central differences, and dw at z0."""
    h = 1e-5 * max(1.0, abs(z0))

    def k(z):
        return candidates.evaluate(key, np.array([z]))[0][0].real

    _, dw = candidates.evaluate(key, np.array([z0]))
    return (k(z0 + h) - k(z0 - h)) / (2 * h), (k(z0 + 1j * h) - k(z0 - 1j * h)) / (2 * h), dw[0]


@pytest.mark.parametrize(("key", "z0"), FORMULA_POINTS)
def test_derivative_is_that_of_the_formula(key, z0):
    # dw = K_x − i K_y of the candidate's own formula, in both directions.
    k_x, k_y, dw = k_x_and_k_y(key, z0)
    assert abs(k_x - dw.real) / abs(dw) < 1e-8
    assert abs(k_y + dw.imag) / abs(dw) < 1e-8
