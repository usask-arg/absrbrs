"""The Faddeeva generator: known values, the zeros at x = 0, an independent check of w′, the
grid of spec §6.2 and the record layout."""

from flint import acb, arb, ctx

from absrbrs_reference.generators import faddeeva


def close(v: arb, expected: float, rel: float = 2.0**-52) -> bool:
    return abs(float(v.mid()) - expected) <= rel * abs(expected)


def test_known_values():
    w, dw, _ = faddeeva.w_and_derivative(0.0, 0.0)
    assert w.real == 1 and w.imag == 0 and dw.real == 0
    assert close(dw.imag, 1.1283791670955126)  # 2/√π
    w, _, _ = faddeeva.w_and_derivative(0.0, 1.0)
    assert close(w.real, 0.42758357615580700)  # erfcx(1)
    w, _, _ = faddeeva.w_and_derivative(1.0, 0.0)
    assert close(w.real, 0.36787944117144233)  # exp(−1)
    assert close(w.imag, 0.6071577058413937)  # (2/√π) D(1), D Dawson's integral


def test_w_is_real_on_the_imaginary_axis():
    w, dw, _ = faddeeva.w_and_derivative(0.0, 2.0)
    assert w.imag == 0 and dw.real == 0


def test_derivative_agrees_with_an_arb_central_difference():
    def w_at(z: acb) -> acb:
        return (-(z * z)).exp() * (-acb(0, 1) * z).erfc()

    for x, y in [(1.5, 0.25), (20.0, 1e-20), (1e6, 1.0)]:
        _, dw, _ = faddeeva.w_and_derivative(x, y)
        with ctx.workprec(1024):
            h = arb(2) ** -150
            z = acb(x, y)
            fd = (w_at(z + h) - w_at(z - h)) / (2 * h)
            assert abs(fd - dw) <= abs(dw) * arb(2) ** -70


def test_grid():
    xs, ys = faddeeva.x_grid(), faddeeva.y_grid()
    assert len(xs) == 124 and xs == sorted(set(xs))
    assert xs[0] == 0.0 and xs[-1] == 1e8 and 1e-8 in xs and 30.0 in xs
    assert len(ys) == 74
    assert ys[0] == 0.0 and ys[1] == 1e-30 and ys[-1] == 1e6


def test_record_layout():
    precision, records = faddeeva.generate(xs=[0.0, 1.0], ys=[0.0, 1.0])
    assert precision >= faddeeva.START_BITS
    assert [r["case"] for r in records] == ["y=0.0", "y=1.0"]
    for r in records:
        assert r["generator"] == "faddeeva" and r["label"] == "certified"
        assert r["inputs"]["x"] == [0.0, 1.0]
        assert set(r["outputs"]) == {"w_re", "w_im", "dw_re", "dw_im"}
        assert all(len(v["hi"]) == 2 for v in r["outputs"].values())
