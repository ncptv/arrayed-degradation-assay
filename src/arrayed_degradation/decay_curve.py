import typing as tp
import warnings

import numpy as np
import numpy.typing as npt
from scipy.optimize import OptimizeWarning, curve_fit
from sklearn.metrics import r2_score

warnings.simplefilter("ignore", OptimizeWarning)


def clip(value: float, min_bound: float = -1e4, max_bound: float = 1e4) -> float:
    return min(max(value, min_bound), max_bound)


def fit_decay_curve(
    durations: npt.NDArray[np.float64], observations: npt.NDArray[np.float64]
) -> dict[str, float | tp.Callable]:
    """
    Fit a decay curve to the given data.

    Args:
        durations (npt.NDArray[np.float64]):
            The time points at which the observations were made.
        observations (npt.NDArray[np.float64]):
            The observed values.

    Returns:
        A dictionary containing the following keys:
        - decay_rate: The estimated decay rate.
        - decay_rate_std: The standard deviation of the estimated decay rate.
        - half_life: The estimated half-life.
        - half_life_std: The standard deviation of the estimated half-life.
        - fit: A function that can be used to evaluate the fitted decay curve.
    """

    popt, pcov = curve_fit(
        lambda x, A, tau: np.exp(-tau * x + A), durations, observations, method="dogbox"
    )
    A_hat = popt[0]
    tau_hat = popt[1]
    fit_fn = lambda x: np.exp(-tau_hat * x + A_hat)
    tau_var = pcov[1, 1]
    half_life = np.log(2) / tau_hat
    half_life_std = np.abs(np.log(2) / tau_hat**2 * np.sqrt(tau_var))

    return {
        "decay_rate": clip(tau_hat),
        "decay_rate_std": clip(np.sqrt(tau_var)),
        "half_life": clip(half_life),
        "half_life_std": clip(half_life_std),
        "fit": fit_fn,
    }


def compute_cv_r2(
    durations: npt.NDArray[np.float64],
    observations: npt.NDArray[np.float64],
    replicate_indices: npt.NDArray[np.float64],
) -> float:
    """Compute the R2 score using cross-validation.

    If there are multiple replicate indices, the function estimates the R2 goodness of fit score
    by running quasi-cross validation. It fits a decay curve to the data excluding one replicate
    at a time, and then predicts the observations for the excluded replicate. The R2 score is
    calculated between the predicted observations and the actual observations for the excluded
    replicate. Finally, the average R2 score is returned.

    If there is only one replicate index, the function returns nan as cross-validation is not possible.

    Args:
        durations (npt.NDArray[np.float64]):
            An array of durations. E.g. [0, 3, 6, 0, 3, 6]
        observations (npt.NDArray[np.float64]):
            An array of observations. E.g. [1, 0.5, 0.25, 1, 0.5, 0.25]
        replicate_indices (npt.NDArray[np.float64]):
            An array of replicate indices. E.g. [0, 0, 0, 1, 1, 1]

    Returns:
        float: The average cross-validated R2 score.

    """
    r2s = []
    # estimate R2 goodness of fit score by running quasi-cross validation
    if len(np.unique(replicate_indices)) > 1:
        for repl_ind in np.unique(replicate_indices):
            mask = replicate_indices != repl_ind
            fit = fit_decay_curve(durations[mask], observations[mask])
            pred = tp.cast(tp.Callable, fit["fit"])(durations[~mask])
            r2s.append(r2_score(observations[~mask], pred))
        return np.array(r2s).mean()
    else:
        return np.nan


def estimate_half_life(
    durations: npt.NDArray[np.float64],
    observations: npt.NDArray[np.float64],
    replicate_indices: npt.NDArray[np.float64],
) -> dict[str, tp.Any]:
    """
    Estimate the half-life of a decay curve.

    Args:
        durations (npt.NDArray[np.float64]):
            An array of durations.
        observations (npt.NDArray[np.float64]):
            An array of observations.
        replicate_indices (npt.NDArray[np.float64]):
            An array of replicate indices.

    Returns:
        dict[str, tp.Any]:
            A dictionary containing the estimated decay rate, decay rate standard deviation,
            half-life, half-life standard deviation, R-squared score, and the fit parameters of the decay curve.
    """

    r2 = compute_cv_r2(durations, observations, replicate_indices)
    fit = fit_decay_curve(durations, observations)
    return {
        "decay_rate": fit["decay_rate"],
        "decay_rate_std": fit["decay_rate_std"],
        "half_life": fit["half_life"],
        "half_life_std": fit["half_life_std"],
        "r2_score": r2,
        "fit": fit["fit"],
    }
