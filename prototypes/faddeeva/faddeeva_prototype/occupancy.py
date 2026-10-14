"""Where a line's window puts its grid points in (|x|, y): the cost weights of the search.

Six representative lines of the primary workload (spec §10.1): an O₂ A-band and a CO₂ 1.6 µm
line, each at the surface, near 30 km and near 65 km, on a 0.005 cm⁻¹ grid with the default
25 cm⁻¹ window. Each line contributes equally.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from faddeeva_prototype.reference import complex_points

K_B = 1.380649e-23  # J/K
N_A = 6.02214076e23  # 1/mol
C = 299792458.0  # m/s

SPACING_CMINV = 0.005
WIDTH_CMINV = 25.0


@dataclass(frozen=True)
class Line:
    """A representative line: centre (cm⁻¹), isotopologue mass (g/mol), γ_air (cm⁻¹/atm), n_air."""

    name: str
    centre: float
    mass: float
    gamma_air: float
    n_air: float


LINES = (
    Line("O2 A-band", 13100.0, 31.98983, 0.045, 0.72),
    Line("CO2 1.6 um", 6300.0, 43.98983, 0.070, 0.75),
)
#: (temperature K, pressure Pa): the surface, about 30 km, about 65 km.
CONDITIONS = ((250.0, 101325.0), (225.0, 1000.0), (200.0, 10.0))


def doppler_width(centre: float, temperature: float, mass: float) -> float:
    """Doppler 1/e half width σ_D = (ν₀/c) √(2 k_B T / m), in cm⁻¹."""
    return centre / C * np.sqrt(2 * K_B * temperature * N_A * 1e3 / mass)


def lorentz_width(line: Line, temperature: float, pressure: float) -> float:
    """Air-broadened half width (296/T)^n γ_air p, in cm⁻¹ (p in atm)."""
    return line.gamma_air * (296.0 / temperature) ** line.n_air * pressure / 101325.0


def workload_points() -> NDArray[np.complex128]:
    """|x| + iy of every grid point in the windows of the representative lines."""
    n = round(WIDTH_CMINV / SPACING_CMINV)
    offsets = np.arange(-n, n + 1) * SPACING_CMINV
    out = []
    for line in LINES:
        for temperature, pressure in CONDITIONS:
            sd = doppler_width(line.centre, temperature, line.mass)
            y = lorentz_width(line, temperature, pressure) / sd
            out.append(complex_points(np.abs(offsets) / sd, np.full(offsets.size, y)))
    return np.concatenate(out)
