"""The Faddeeva function w(z) = exp(−z²) erfc(−iz) and its derivative, on the validation grid.

One record per y value (case ``y=<repr(y)>``): inputs ``x`` (doubles ≥ 0, increasing) and ``y``
(a double ≥ 0); outputs ``w_re``, ``w_im``, ``dw_re`` and ``dw_im`` at each x. Negative x follow
from the symmetries w(−x + iy) = conj(w(x + iy)) and w′(−x + iy) = −conj(w′(x + iy)) (so Re w
and Im w′ are even in x, Im w and Re w′ odd) and are not stored.

w′(z) = −2z w(z) + 2i/√π is formed in Arb, where the cancellation this identity suffers in
floating point (about |z|² relative) only costs working precision. Each point is recomputed
at doubling precision until every stored component is known to 2^-80 (``TARGET_BITS``) relative
to the scales of the spec §6.2 metrics: Re w relative to itself, Im w relative to |w|, and both
parts of w′ relative to |w′|. Values below the double range (Re w at y = 0, x ≥ 27) are stored
as subnormals or 0, with an absolute ``err``.
At x = 0, w(iy) is real, so Im w and Re w′ are exactly zero; they are stored as exact zeros.

The grid (spec §6.2): x ∈ {0, 0.5, …, 30} or {10^(k/4) : |k| ≤ 32}, 1e-8 to 1e8; y ∈ {0} or
{10^(k/2) : −60 ≤ k ≤ 12}, 1e-30 to 1e6.
"""

from __future__ import annotations

from flint import acb, arb, ctx

from absrbrs_reference.records import nearest_pow10, record

TARGET_BITS = 80
START_BITS = 128
MAX_BITS = 16384


def x_grid() -> list[float]:
    """0, 0.5, …, 30 (exact) and 10^(k/4) for |k| ≤ 32, sorted and without duplicates."""
    linear = [k / 2 for k in range(61)]
    logs = [nearest_pow10(k / 4) for k in range(-32, 33)]
    return sorted(set(linear) | set(logs))


def y_grid() -> list[float]:
    """0 and 10^(k/2) for −60 ≤ k ≤ 12."""
    return [0.0, *(nearest_pow10(k / 2) for k in range(-60, 13))]


def _accurate(part: arb, scale: arb) -> bool:
    """Whether the ball ``part`` is known to TARGET_BITS relative to ``scale`` (> 0)."""
    return bool(part.rad() <= scale * arb(2) ** -TARGET_BITS)


def w_and_derivative(x: float, y: float) -> tuple[acb, acb, int]:
    """Arb balls for w(x + iy) and w′(x + iy), and the precision (bits) that achieved them.

    Re w is known to TARGET_BITS relative to itself (it is positive), Im w relative to |w|, and
    both parts of w′ relative to |w′|, the scales of the spec §6.2 metrics. Raises
    ``ArithmeticError`` if MAX_BITS is not enough.
    """
    bits = START_BITS
    while bits <= MAX_BITS:
        with ctx.workprec(bits):
            z = acb(x, y)
            w = (-(z * z)).exp() * (-acb(0, 1) * z).erfc()
            dw = -2 * z * w + 2 * acb(0, 1) / arb.pi().sqrt()
            if x == 0.0:
                w = acb(w.real, 0)
                dw = acb(0, dw.imag)
            w_scale, dw_scale = abs(w).mid(), abs(dw).mid()
            if (
                _accurate(w.real, abs(w.real.mid()))
                and (x == 0.0 or _accurate(w.imag, w_scale))
                and (x == 0.0 or _accurate(dw.real, dw_scale))
                and _accurate(dw.imag, dw_scale)
            ):
                return w, dw, bits
        bits *= 2
    raise ArithmeticError(f"w({x!r} + {y!r}i) not resolved to 2^-{TARGET_BITS} by {MAX_BITS} bits")


def generate(
    xs: list[float] | None = None, ys: list[float] | None = None
) -> tuple[int, list[dict]]:
    """Records for every y of ``ys`` (default `y_grid`) over the x of ``xs`` (default `x_grid`)."""
    xs = x_grid() if xs is None else xs
    ys = y_grid() if ys is None else ys
    records = []
    max_bits = START_BITS
    for y in ys:
        outputs: dict[str, list[arb]] = {"w_re": [], "w_im": [], "dw_re": [], "dw_im": []}
        for x in xs:
            w, dw, bits = w_and_derivative(x, y)
            max_bits = max(max_bits, bits)
            outputs["w_re"].append(w.real)
            outputs["w_im"].append(w.imag)
            outputs["dw_re"].append(dw.real)
            outputs["dw_im"].append(dw.imag)
        records.append(record("faddeeva", f"y={y!r}", {"x": xs, "y": y}, outputs))
    return max_bits, records
