"""Candidate approximations of the Faddeeva function w(z) = exp(-z²) erfc(-iz), for Im z ≥ 0.

Every candidate returns ``(w, dw)``: w(z) and its derivative w′(z), both from the same formula
(the derivative of the approximation itself, never the identity w′ = −2zw + 2i/√π), as
complex128 arrays shaped like ``z``.
"""

from __future__ import annotations

from functools import cache
from math import factorial

import numpy as np
from numpy.typing import NDArray
from scipy.special import dawsn

Complex = NDArray[np.complex128]

INV_SQRT_PI = 0.5641895835477562869480795

# sasktran2's 8-pole J-pole constants (rust/sasktran2-rs/src/math/errorfunctions/optimized.rs).
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
    return a[1 : n + 1][::-1].copy(), float(big_l)


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
    depth 1 is Humlíček's region-I formula iz/(√π (z² − ½))."""
    f = np.array(z, dtype=np.complex128, copy=True)
    fp = np.ones_like(f)
    for k in range(depth, 0, -1):
        a = k / 2
        f, fp = z - a / f, 1 + a * fp / f**2
    return 1j * INV_SQRT_PI / f, -1j * INV_SQRT_PI * fp / f**2


def gauss_cf(z: Complex, depth: int) -> tuple[Complex, Complex]:
    """exp(−z²) plus the continued fraction: keeps the Gaussian core that the asymptotic
    expansion lacks, for small y and moderate |x|."""
    g = np.exp(-(z * z))
    w, dw = laplace_cf(z, depth)
    return w + g, dw - 2 * z * g


def small_y_taylor(z: Complex, order: int) -> tuple[Complex, Complex]:
    """Taylor series in iy about the real axis, to (iy)^order:
    w(x + iy) ≈ Σₖ (iy)ᵏ w⁽ᵏ⁾(x)/k!, with w(x) = exp(−x²) + (2i/√π) D(x) (D: Dawson) and
    w⁽ᵏ⁺¹⁾ = −2x w⁽ᵏ⁾ − 2k w⁽ᵏ⁻¹⁾."""
    x = z.real
    iy = 1j * z.imag
    w0 = np.exp(-(x * x)) + 2j * INV_SQRT_PI * dawsn(x)
    derivs = [w0, -2 * x * w0 + 2j * INV_SQRT_PI]
    for k in range(1, order + 1):
        derivs.append(-2 * x * derivs[k] - 2 * k * derivs[k - 1])
    w = np.zeros_like(w0)
    dw = np.zeros_like(w0)
    for k in range(order + 1):
        c = iy**k / factorial(k)
        w = w + c * derivs[k]
        dw = dw + c * derivs[k + 1]
    return w, dw


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
