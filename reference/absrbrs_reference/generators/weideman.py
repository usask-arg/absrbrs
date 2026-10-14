"""Weideman's (1994) N = 24 rational approximation of w(z): the constant L and coefficients aₙ.

w(z) ≈ 2 Σₙ aₙ Zⁿ⁻¹/(L − iz)² + 1/(√π (L − iz)), Z = (L + iz)/(L − iz), L = √(N/√2).

The aₙ are the discrete Fourier coefficients of F(t) = exp(−t²)(L² + t²) at the 2M = 4N points
t_k = L tan(kπ/(2M)): aₙ = (F(0) + 2 Σₖ₌₁^(M−1) F(t_k) cos(πkn/M)) / (2M), n = 1 … N. (F(t_M) is
taken as 0, the limit as t → ∞.) This is the formula of the paper; its code listing is not used.
"""

from __future__ import annotations

from flint import arb, ctx

from absrbrs_reference.records import record

N = 24
PRECISION_BITS = 256


def coefficients() -> tuple[arb, list[arb]]:
    """L and [a₁, …, a_N] as Arb balls."""
    m = 2 * N
    big_l = (arb(N) / arb(2).sqrt()).sqrt()

    def f(k: int) -> arb:
        t = big_l * (arb.pi() * k / (2 * m)).tan()
        return (-(t * t)).exp() * (big_l * big_l + t * t)

    values = [f(k) for k in range(m)]
    a = [
        (values[0] + 2 * sum(values[k] * (arb.pi() * k * n / m).cos() for k in range(1, m)))
        / (2 * m)
        for n in range(1, N + 1)
    ]
    return big_l, a


def generate() -> tuple[int, list[dict]]:
    with ctx.workprec(PRECISION_BITS):
        big_l, a = coefficients()
        return PRECISION_BITS, [
            record("weideman", "n_terms=24", {"n_terms": N}, {"L": [big_l], "a": a})
        ]
