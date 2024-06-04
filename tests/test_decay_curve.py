import numpy as np
import pytest
from pytest import approx

from arrayed_degradation.decay_curve import (
    estimate_half_life,
)

excellent_decay_data = (
    np.array([0, 3, 6, 0, 3, 6, 0, 3, 6]),
    np.array([1, 0.5, 0.25, 1, 0.5, 0.25, 1, 0.5, 0.25]),
    np.array([0, 0, 0, 1, 1, 1, 2, 2, 2]),
    {
        "decay_rate": 0.231049,
        "decay_rate_std": 0,
        "half_life": 3,
        "half_life_std": 0,
        "fit": lambda x: np.exp(-0.231049 * x),
    },
)

noisy_decay_data = (
    np.array([0, 3, 6, 0, 3, 6, 0, 3, 6]),
    np.array([1, 0.51, 0.24, 1, 0.48, 0.24, 1, 0.52, 0.26]),
    np.array([0, 0, 0, 1, 1, 1, 2, 2, 2]),
    {
        "decay_rate": 0.231250184,
        "decay_rate_std": 0.00401658,
        "half_life": 2.99739,
        "half_life_std": 0.0520616,
        "fit": lambda x: np.exp(-0.231250184 * x + 0.00080667),
    },
)


@pytest.mark.parametrize(
    "durations,observations,replicate_indices,expected_result",
    [excellent_decay_data, noisy_decay_data],
)
def test_estimate_half_life(
    durations, observations, replicate_indices, expected_result
):
    result = estimate_half_life(durations, observations, replicate_indices)

    assert result["half_life"] == approx(expected_result["half_life"])
    assert result["half_life_std"] == approx(expected_result["half_life_std"])
    assert result["decay_rate_std"] == approx(expected_result["decay_rate_std"])
    assert result["decay_rate"] == approx(expected_result["decay_rate"])
    assert result["fit"](0) == approx(expected_result["fit"](0))
    assert result["fit"](3) == approx(expected_result["fit"](3))
