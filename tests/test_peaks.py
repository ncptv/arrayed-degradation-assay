import numpy as np
import pytest
from pytest import approx

from arrayed_degradation.input_data import EPG
from arrayed_degradation.peaks import (
    compute_peak_bounds,
    get_peak_area,
    normalize_trace_wrt_control,
)
from arrayed_degradation.tests.utils import gaussian_pdf


@pytest.fixture
def clean_epg() -> EPG:
    x = np.linspace(0, 1000, 1001)
    trace = gaussian_pdf(x, 500, 50)
    trace = trace / np.sum(trace)
    return EPG(nucleotides=x, trace=trace)


@pytest.fixture
def epg_with_background() -> EPG:
    # after background subtraction the peak area should be 1
    x = np.linspace(0, 1000, 1001)
    trace = gaussian_pdf(x, 500, 50)
    trace = trace / np.sum(trace) + 2
    return EPG(nucleotides=x, trace=trace)


@pytest.fixture
def noisy_epg() -> EPG:
    x = np.linspace(0, 1000, 1001)
    y_target = gaussian_pdf(x, 500, 50)
    y_secondary = gaussian_pdf(x, 360, 15) * 0.1
    return EPG(nucleotides=x, trace=y_target + y_secondary)


@pytest.fixture
def epg_with_control_peak() -> EPG:
    x = np.linspace(0, 1000, 1001)
    y_target = gaussian_pdf(x, 500, 50) * 3
    y_control = gaussian_pdf(x, 150, 6) * 0.5
    return EPG(nucleotides=x, trace=y_target + y_control)


@pytest.mark.parametrize(
    "epg, rel_height, expected_left, expected_right",
    [
        ("clean_epg", 0.75, 416, 583),
        ("clean_epg", 0.9, 392, 607),
        ("noisy_epg", 0.75, 416, 583),
        ("noisy_epg", 0.9, 336, 607),
    ],
)
def test_peak_bounds(
    request, epg, rel_height, expected_left: float, expected_right: float
):
    min_peak = 0
    max_peak = 1000
    epg = request.getfixturevalue(epg)
    left, peak, right, prominence = compute_peak_bounds(
        epg.nucleotides, epg.trace, min_peak, max_peak, rel_height
    )
    assert left == approx(expected_left)
    assert right == approx(expected_right)
    assert peak == approx(500)


@pytest.mark.parametrize(
    "epg, min_peak, max_peak, remove_background, prominence, expected_area",
    [
        ("clean_epg", 416, 583, False, None, 0.905053),
        ("clean_epg", 416, 583, True, None, 0.574612),
        ("clean_epg", 0, 1000, False, None, 1),
        ("epg_with_background", 0, 1000, True, 0.0079788, 1),
        ("noisy_epg", 416, 583, False, None, 0.9050626),
        ("noisy_epg", 416, 583, True, None, 0.574413),
    ],
)
def test_get_peak_area(
    request,
    epg,
    min_peak: float,
    max_peak: float,
    remove_background: bool,
    prominence: float | None,
    expected_area: float,
):
    epg = request.getfixturevalue(epg)
    area = get_peak_area(
        epg.nucleotides, epg.trace, min_peak, max_peak, remove_background, prominence
    )
    assert area == approx(expected_area, abs=1e-3)


def test_normalize_trace_wrt_control(epg_with_control_peak: EPG):
    control_peak_min = 100
    control_peak_max = 200
    normalized_epg = normalize_trace_wrt_control(
        epg_with_control_peak, control_peak_min, control_peak_max
    )
    assert normalized_epg.trace[control_peak_min:control_peak_max].sum() == approx(
        1, abs=1e-3
    )
    assert (epg_with_control_peak.trace / normalized_epg.trace) == approx(0.5)
    assert normalized_epg.trace.sum() == approx(7)
