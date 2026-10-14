"""The line-state quantities q = (A, ν₀′, σ_D, y) of one line at one condition (spec §5.2).

Pressures are in atm: p_atm = p/101325, p_self = VMR·p_atm, p_air = (1 − VMR)·p_atm.

* f = S(T)/S_ref = exp(L) · exp(−c₂E″(1/T − 1/T_ref)) · expm1(−c₂ν₀/T)/expm1(−c₂ν₀/T_ref), with
  L = ln(Q(T_ref)/Q(T)) given as an input;
* σ_D = (ν₀/c) √(2k_B T/m), m = mass·10⁻³/N_A kg;
* γ_L = (T_ref/T)^{n_air} γ_air p_air + (T_ref/T)^{n_self} γ_self p_self;
* ν₀′ = ν₀ + (δ_air + δ′_air (T − T_ref)) p_air + δ_self p_self;
* A = 10⁻⁴ S_ref f/(√π σ_D) (m²/molecule) and y = γ_L/σ_D.

An absent n_self (``null``) is n_air and an absent δ_self is δ_air. The constants are the exact
SI values, so c₂ = hc/k_B is an exact rational; every input is a double, converted to Arb
exactly. One record per (line, T, p, VMR): case ``line=<k>,T=<T>,p=<p>,vmr=<VMR>``.
"""

from __future__ import annotations

from fractions import Fraction

from flint import arb, ctx, fmpq

from absrbrs_reference.records import record

PRECISION_BITS = 256

# Exact SI constants.
H = Fraction(662607015, 10**42)  # J s
C = Fraction(299792458)  # m/s
K_B = Fraction(1380649, 10**29)  # J/K
N_A = Fraction(602214076) * 10**15  # 1/mol
C2 = H * C / K_B * 100  # cm K
T_REF = Fraction(296)
ATM = Fraction(101325)

#: (ν₀, S_ref, E″, γ_air, γ_self, n_air, n_self, δ_air, δ′_air, δ_self); None is absent.
LINES = [
    (13100.5, 7.4e-24, 1200.0, 0.045, 0.047, 0.72, 0.70, -0.008, 1e-5, -0.007),
    (6300.25, 1.6e-23, 300.0, 0.07, 0.09, 0.75, None, -0.004, 0.0, None),
    (0.75, 1e-24, 45.0, 0.1, 0.4, 0.68, None, 0.0, 0.0, None),  # microwave: c₂ν₀/T small
    (2140.0, 4e-19, 8000.0, 0.05, 0.06, 0.7, 0.6, -0.002, 0.0, 0.003),
]
FIELDS = (
    "nu",
    "s_ref",
    "e_lower",
    "gamma_air",
    "gamma_self",
    "n_air",
    "n_self",
    "delta_air",
    "deltap_air",
    "delta_self",
)
TEMPERATURES = [50.0, 150.0, 296.0, 1500.0]
PRESSURES = [0.0, 1000.0, 101325.0, 1e6]
VMRS = [0.0, 0.01, 1.0]
LOG_Q_RATIOS = [0.0, 0.37, -1.2]  # one per line, cycling
MASSES = [32.0, 44.0, 18.0, 28.0]  # g/mol, one per line


def exact(x: float | Fraction) -> arb:
    """``x`` as an exact Arb value (a double or a rational)."""
    q = Fraction(x)
    return arb(fmpq(q.numerator, q.denominator))


def line_shape(
    temperature: float,
    pressure: float,
    vmr: float,
    line: dict[str, float | None],
    log_q_ratio: float,
    mass: float,
) -> dict[str, arb]:
    """q as Arb balls at the current working precision."""
    t, p, x = exact(temperature), exact(pressure), exact(vmr)
    nu, s_ref, e_lower = exact(line["nu"]), exact(line["s_ref"]), exact(line["e_lower"])
    gamma_air, gamma_self = exact(line["gamma_air"]), exact(line["gamma_self"])
    n_air, delta_air, deltap_air = (
        exact(line["n_air"]),
        exact(line["delta_air"]),
        exact(line["deltap_air"]),
    )
    n_self = n_air if line["n_self"] is None else exact(line["n_self"])
    delta_self = delta_air if line["delta_self"] is None else exact(line["delta_self"])
    c2, t_ref = exact(C2), exact(T_REF)

    p_atm = p / exact(ATM)
    p_self = x * p_atm
    p_air = (1 - x) * p_atm
    stim = (-(c2 * nu / t)).expm1() / (-(c2 * nu / t_ref)).expm1()
    boltz = (-(c2 * e_lower) * (1 / t - 1 / t_ref)).exp()
    intensity = exact(log_q_ratio).exp() * boltz * stim
    m = exact(Fraction(mass) / 1000 / N_A)  # kg
    sigma_d = nu / exact(C) * (2 * exact(K_B) * t / m).sqrt()
    ratio = t_ref / t
    gamma = ratio**n_air * gamma_air * p_air + ratio**n_self * gamma_self * p_self
    return {
        "amplitude": exact(Fraction(1, 10000)) * s_ref * intensity / (arb.pi().sqrt() * sigma_d),
        "centre": nu + (delta_air + deltap_air * (t - t_ref)) * p_air + delta_self * p_self,
        "doppler_width": sigma_d,
        "y": gamma / sigma_d,
    }


def generate() -> tuple[int, list[dict]]:
    records = []
    with ctx.workprec(PRECISION_BITS):
        for k, values in enumerate(LINES):
            line = dict(zip(FIELDS, values, strict=True))
            log_q_ratio = LOG_Q_RATIOS[k % len(LOG_Q_RATIOS)]
            mass = MASSES[k]
            for t in TEMPERATURES:
                for p in PRESSURES:
                    for vmr in VMRS:
                        q = line_shape(t, p, vmr, line, log_q_ratio, mass)
                        inputs = {
                            "temperature": t,
                            "pressure": p,
                            "vmr": vmr,
                            **line,
                            "log_q_ratio": log_q_ratio,
                            "mass": mass,
                        }
                        records.append(
                            record(
                                "line_state",
                                f"line={k},T={t!r},p={p!r},vmr={vmr!r}",
                                inputs,
                                {name: [v] for name, v in q.items()},
                                label="certified",
                            )
                        )
    return PRECISION_BITS, records
