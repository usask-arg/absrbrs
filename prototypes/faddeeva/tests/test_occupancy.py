import numpy as np
import pytest

from faddeeva_prototype import occupancy


def test_doppler_width_of_an_o2_a_band_line():
    assert occupancy.doppler_width(13100.0, 250.0, 31.98983) == pytest.approx(0.015752, rel=1e-4)


def test_workload_points_cover_every_window_point_of_every_line():
    z = occupancy.workload_points()
    assert z.size == len(occupancy.LINES) * len(occupancy.CONDITIONS) * 10001
    assert np.all(z.real >= 0) and np.all(z.imag > 0)
