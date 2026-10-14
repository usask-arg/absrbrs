import numpy as np
import pytest

import absrbrs
from absrbrs import _absrbrs, special

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


def test_invalid_point_message_names_the_index_and_value():
    with pytest.raises(absrbrs.InvalidInputError, match=r"z\[1\] = 1-0\.001i: Im z"):
        special.faddeeva([0.0, 1.0 - 1e-3j, 2.0 - 1j])


@pytest.mark.parametrize(
    ("x", "y"),
    [(1.0, -1.0), (np.nan, 1.0), (np.inf, 1.0), (-np.inf, 1.0), (1.0, np.nan), (1.0, np.inf)],
    ids=["negative-y", "nan-x", "inf-x", "-inf-x", "nan-y", "inf-y"],
)
def test_invalid_kernel_points_raise(x, y):
    with pytest.raises(absrbrs.InvalidInputError, match=r"x\[0\] = .*y\[0\] = .*y must be >= 0"):
        special.voigt_kernel(x, y)


def test_negative_zero_is_accepted():
    for tier in ("fast", "accurate"):
        w0, dw0 = special.faddeeva([complex(0.0, 0.0), 2.0 + 0j], tier=tier)
        w1, dw1 = special.faddeeva([complex(-0.0, -0.0), 2.0 - 0j], tier=tier)
        np.testing.assert_array_equal(w1, w0)
        np.testing.assert_array_equal(dw1, dw0)
        for got, expected in zip(
            special.voigt_kernel([-0.0, 2.0], -0.0, tier=tier),
            special.voigt_kernel([0.0, 2.0], 0.0, tier=tier),
            strict=True,
        ):
            np.testing.assert_array_equal(got, expected)


def test_empty_input_gives_empty_outputs():
    w, dw = special.faddeeva([])
    assert w.shape == dw.shape == (0,)
    assert w.dtype == dw.dtype == np.complex128
    k, k_x, k_y = special.voigt_kernel(np.zeros((0, 3)), 1.0)
    assert k.shape == k_x.shape == k_y.shape == (0, 3)
    assert k.dtype == np.float64


def test_memory_layout_does_not_matter():
    rng = np.random.default_rng(2)
    z = rng.uniform(-5, 5, (6, 8)) + 1j * rng.uniform(0, 5, (6, 8))
    x, y = rng.uniform(-5, 5, (6, 8)), rng.uniform(0, 5, (6, 8))
    for view_z, view_x, view_y in [
        (np.asfortranarray(z), np.asfortranarray(x), np.asfortranarray(y)),
        (z[::-1, ::-1], x[::-1, ::-1], y[::-1, ::-1]),
        (z[:, ::-2], x[:, ::-2], y[:, ::-2]),
    ]:
        contiguous_z = np.ascontiguousarray(view_z)
        for got, expected in zip(
            special.faddeeva(view_z), special.faddeeva(contiguous_z), strict=True
        ):
            np.testing.assert_array_equal(got, expected)
        contiguous = (np.ascontiguousarray(view_x), np.ascontiguousarray(view_y))
        for got, expected in zip(
            special.voigt_kernel(view_x, view_y), special.voigt_kernel(*contiguous), strict=True
        ):
            np.testing.assert_array_equal(got, expected)


def test_extension_rejects_unequal_lengths():
    with pytest.raises(absrbrs.InvalidInputError, match="2 and 3 values"):
        _absrbrs.voigt_kernel(np.zeros(2), np.zeros(3), "fast")


def test_unknown_tier_raises():
    with pytest.raises(absrbrs.InvalidInputError, match="tier"):
        special.faddeeva(1.0, tier="exact")


def test_non_numeric_input_raises_type_error():
    with pytest.raises(TypeError):
        special.faddeeva("1.0")
    with pytest.raises(TypeError):
        special.voigt_kernel("1.0", 1.0)


def test_complex_and_bool_arguments_raise_type_error():
    with pytest.raises(TypeError, match="x must be real numbers, not an array of dtype complex128"):
        special.voigt_kernel(1.0 + 0j, 1.0)
    with pytest.raises(TypeError, match="y must be real numbers, not an array of dtype complex128"):
        special.voigt_kernel(1.0, [1.0 + 0j])
    with pytest.raises(TypeError, match="x must be real numbers, not an array of dtype bool"):
        special.voigt_kernel(True, 1.0)
    with pytest.raises(TypeError, match="y must be real numbers, not an array of dtype bool"):
        special.voigt_kernel(1.0, np.array([True]))
    with pytest.raises(TypeError, match="z must be numbers, not an array of dtype bool"):
        special.faddeeva(True)


def test_integer_arguments_are_accepted():
    k, _, _ = special.voigt_kernel([0, 1], 1)
    np.testing.assert_array_equal(k, special.voigt_kernel([0.0, 1.0], 1.0)[0])
    w, _ = special.faddeeva([1, 2])
    np.testing.assert_array_equal(w, special.faddeeva([1.0 + 0j, 2.0 + 0j])[0])
