import numpy as np
import pytest

import absrbrs
from absrbrs import special

SQRT_PI = np.sqrt(np.pi)


def test_known_values():
    for tier in ("fast", "accurate"):
        w, dw = special.faddeeva([0.0, 1.0j, 1.0], tier=tier)
        np.testing.assert_allclose(w.real, [1.0, 0.42758357615580700, np.exp(-1.0)], rtol=1e-6)
        np.testing.assert_allclose(w.imag, [0.0, 0.0, 0.60715770584139372], rtol=1e-6, atol=1e-15)
        np.testing.assert_allclose(dw[0], 2j / SQRT_PI, rtol=1e-6)


def test_tiers_agree_within_the_fast_targets():
    rng = np.random.default_rng(1)
    z = rng.uniform(0, 30, 2000) + 1j * 10.0 ** rng.uniform(-6, 2, 2000)
    w_f, dw_f = special.faddeeva(z, tier="fast")
    w_a, dw_a = special.faddeeva(z, tier="accurate")
    assert np.max(np.abs(w_f.real - w_a.real) / w_a.real) <= 1e-6
    assert np.max(np.abs(dw_f - dw_a) / np.abs(dw_a)) <= 1e-5


def test_fast_is_the_default():
    z = np.array([0.3 + 0.2j, 12.0 + 1e-3j])
    for got, expected in zip(special.faddeeva(z), special.faddeeva(z, tier="fast"), strict=True):
        np.testing.assert_array_equal(got, expected)


def test_shapes_broadcast_and_scalars_give_zero_d_arrays():
    w, dw = special.faddeeva(np.ones((2, 3)) * (1 + 1j))
    assert w.shape == dw.shape == (2, 3)
    w, _ = special.faddeeva(1.0)
    assert w.shape == ()
    k, k_x, k_y = special.voigt_kernel(np.linspace(0, 5, 4)[:, None], [0.0, 0.5, 2.0])
    assert k.shape == k_x.shape == k_y.shape == (4, 3)


def test_voigt_kernel_is_the_triplet_of_w():
    x, y = np.linspace(-20, 20, 81), 0.25
    k, k_x, k_y = special.voigt_kernel(x, y)
    w, dw = special.faddeeva(x + 1j * y)
    np.testing.assert_array_equal(k, w.real)
    np.testing.assert_array_equal(k_x, dw.real)
    np.testing.assert_array_equal(k_y, -dw.imag)


def test_profile_normalization():
    # ∫ K dx = √π for any y.
    x = np.linspace(-2000, 2000, 400_001)
    k, _, _ = special.voigt_kernel(x, 1.0)
    integral = (
        np.sum((k[1:] + k[:-1]) * np.diff(x)) / 2
    )  # trapezoid rule (numpy 1.26 has no trapezoid)
    assert integral == pytest.approx(SQRT_PI, rel=1e-3)


@pytest.mark.parametrize(
    "z", [1.0 - 1e-3j, complex(1.0, np.nan), complex(np.inf, 1.0)], ids=["negative", "nan", "inf"]
)
def test_invalid_points_raise(z):
    with pytest.raises(absrbrs.InvalidInputError, match="Im z"):
        special.faddeeva([0.0, z])


def test_invalid_kernel_points_raise():
    with pytest.raises(absrbrs.InvalidInputError):
        special.voigt_kernel(1.0, -1.0)


def test_unknown_tier_raises():
    with pytest.raises(absrbrs.InvalidInputError, match="tier"):
        special.faddeeva(1.0, tier="exact")


def test_non_numeric_input_raises_type_error():
    with pytest.raises(TypeError):
        special.faddeeva("1.0")
