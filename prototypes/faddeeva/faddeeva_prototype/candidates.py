"""Candidate approximations of the Faddeeva function w(z) = exp(-z²) erfc(-iz), for Im z ≥ 0.

Every candidate returns ``(w, dw)`` as complex128 arrays shaped like ``z``: w(z), and
dw = K_x − i K_y, where K_x and K_y are the x- and y-derivatives of K = Re w *of the candidate's
own formula*. For the analytic candidates dw is the derivative w′(z) of the approximation (never
the identity w′ = −2zw + 2i/√π applied to it); for the Taylor series, which is not analytic in
z, it is built from the formula's two partial derivatives.

Each candidate is meant for part of the upper half-plane only, as its docstring says; the
composites of `search` use them only there.
"""

from __future__ import annotations

from functools import cache
from math import factorial

import numpy as np
from numpy.typing import NDArray
from scipy.special import dawsn

from faddeeva_prototype.reference import complex_points

Complex = NDArray[np.complex128]

INV_SQRT_PI = 0.5641895835477562869480795

# sasktran2's 8-pole J-pole constants. From sasktran2
# (rust/sasktran2-rs/src/math/errorfunctions/optimized.rs): Copyright (c) 2024 University of
# Saskatchewan Atmospheric Research Group, MIT licence.
_JPOLE_B = np.array(
    [
        0.00383968430671409 - 0.0119854387180615j,
        -0.321597857664957 - 0.218883985607935j,
        2.55515264319988 + 0.613958600684469j,
        -2.73739446984183 + 5.69007914897806j,
        0.00383968430671409 + 0.0119854387180615j,
        -0.321597857664957 + 0.218883985607935j,
        2.55515264319988 - 0.613958600684469j,
        -2.73739446984183 - 5.69007914897806j,
    ]
)
_JPOLE_C = np.array(
    [
        2.51506776338386 - 1.60713668042405j,
        -1.68985621846204 - 1.66471695485661j,
        0.981465428659098 - 1.70017951305004j,
        -0.322078795578047 - 1.71891780447016j,
        -2.51506776338386 - 1.60713668042405j,
        1.68985621846204 - 1.66471695485661j,
        -0.981465428659098 - 1.70017951305004j,
        0.322078795578047 - 1.71891780447016j,
    ]
)


def jpole8(z: Complex) -> tuple[Complex, Complex]:
    """sasktran2's 8-pole approximation (the baseline): w ≈ −(i/√π) Σⱼ bⱼ/(z − cⱼ)."""
    d = z[..., None] - _JPOLE_C
    w = (-1j * INV_SQRT_PI) * np.sum(_JPOLE_B / d, axis=-1)
    dw = (1j * INV_SQRT_PI) * np.sum(_JPOLE_B / d**2, axis=-1)
    return w, dw


@cache
def weideman_coefficients(n: int) -> tuple[NDArray[np.float64], float]:
    """Weideman (1994) coefficients ``a`` (highest power first, for np.polyval) and L."""
    m = 2 * n
    k = np.arange(-m + 1, m)
    big_l = np.sqrt(n / np.sqrt(2.0))
    t = big_l * np.tan(k * np.pi / (2 * m))
    f = np.concatenate([[0.0], np.exp(-(t**2)) * (big_l**2 + t**2)])
    a = np.real(np.fft.fft(np.fft.fftshift(f))) / (2 * m)
    a = a[1 : n + 1][::-1].copy()
    a.flags.writeable = False  # shared by every caller through the cache
    return a, float(big_l)


def weideman(z: Complex, n: int) -> tuple[Complex, Complex]:
    """Weideman's rational approximation with n terms:
    w = 2 p(Z)/(L − iz)² + 1/(√π (L − iz)), Z = (L + iz)/(L − iz)."""
    a, big_l = weideman_coefficients(n)
    d = big_l - 1j * z
    zz = (big_l + 1j * z) / d
    p = np.polyval(a, zz)
    dp = np.polyval(np.polyder(a), zz)
    w = 2 * p / d**2 + INV_SQRT_PI / d
    dw = 2 * dp * (2j * big_l / d**2) / d**2 + 4j * p / d**3 + 1j * INV_SQRT_PI / d**2
    return w, dw


def laplace_cf(z: Complex, depth: int) -> tuple[Complex, Complex]:
    """Laplace continued fraction truncated after ``depth`` levels:
    w ≈ (i/√π) / (z − ½/(z − 1/(z − (3/2)/(z − …)))). Depth 0 is the Lorentzian i/(√π z);
    depth 1 is Humlíček's region-I formula iz/(√π (z² − ½)).

    An asymptotic expansion: accurate for large |z| only. The truncated fraction has poles on
    the real axis (NaN there), all within |x| < 6 for depth ≤ 16; use it for |z| ≥ 8."""
    f = np.array(z, dtype=np.complex128, copy=True)
    fp = np.ones_like(f)
    for k in range(depth, 0, -1):
        a = k / 2
        f, fp = z - a / f, 1 + a * fp / f**2
    return 1j * INV_SQRT_PI / f, -1j * INV_SQRT_PI * fp / f**2


def gauss_cf(z: Complex, depth: int) -> tuple[Complex, Complex]:
    """exp(−z²) plus the continued fraction: keeps the Gaussian core that the asymptotic
    expansion lacks, for small y and moderate |x|. exp(−z²) = exp(y² − x²) e^(−2ixy) overflows
    for y above about 26.6, so this is for small y only."""
    g = np.exp(-(z * z))
    w, dw = laplace_cf(z, depth)
    return w + g, dw - 2 * z * g


def small_y_taylor(z: Complex, order: int) -> tuple[Complex, Complex]:
    """Taylor series in iy about the real axis, to (iy)^order:
    w(x + iy) ≈ Σₖ (iy)ᵏ w⁽ᵏ⁾(x)/k!, with w(x) = exp(−x²) + (2i/√π) D(x) (D: Dawson) and
    w⁽ᵏ⁺¹⁾ = −2x w⁽ᵏ⁾ − 2k w⁽ᵏ⁻¹⁾.

    The truncated series is not analytic in z, so its y-derivative is not i times its
    x-derivative: dw is built from ∂K/∂x = Re Σₖ (iy)ᵏ w⁽ᵏ⁺¹⁾/k! and
    ∂K/∂y = Re(i Σₖ₍ₖ₌₀…ₙ₋₁₎ (iy)ᵏ w⁽ᵏ⁺¹⁾/k!) of the formula itself. On the real axis w′(x)
    comes from the identity −2xw + 2i/√π, which loses about 2x² in relative accuracy: for
    |x| ≲ 6 only."""
    x = z.real
    iy = 1j * z.imag
    w0 = np.exp(-(x * x)) + 2j * INV_SQRT_PI * dawsn(x)
    derivs = [w0, -2 * x * w0 + 2j * INV_SQRT_PI]
    for k in range(1, order + 1):
        derivs.append(-2 * x * derivs[k] - 2 * k * derivs[k - 1])
    w = np.zeros_like(w0)
    d_x = np.zeros_like(w0)
    d_y = np.zeros_like(w0)
    for k in range(order + 1):
        c = iy**k / factorial(k)
        w = w + c * derivs[k]
        d_x = d_x + c * derivs[k + 1]
        if k < order:
            d_y = d_y + c * derivs[k + 1]
    k_x = d_x.real
    k_y = (1j * d_y).real
    return w, complex_points(k_x, -k_y)


#: Every method the search may use, by key ``(name, parameter)``.
METHODS = {
    "jpole8": lambda z, _: jpole8(z),
    "weideman": weideman,
    "cf": laplace_cf,
    "gauss_cf": gauss_cf,
    "taylor": small_y_taylor,
}


def evaluate(key: tuple[str, int], z: Complex) -> tuple[Complex, Complex]:
    """Evaluate the method ``key`` (for example ``("weideman", 32)``) at ``z``."""
    name, param = key
    return METHODS[name](z, param)
