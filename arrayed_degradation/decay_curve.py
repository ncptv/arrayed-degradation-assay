import typing as tp
import warnings

import numpy as np
import numpy.typing as npt
from scipy.optimize import curve_fit, OptimizeWarning
from sklearn.metrics import r2_score

warnings.simplefilter("ignore", OptimizeWarning)


def clip(value: float, min_bound: float = -1e6, max_bound: float = 1e6) -> float:
    return min(max(value, min_bound), max_bound)


def fit_decay_curve(
    durations: npt.NDArray[np.float64], observations: npt.NDArray[np.float64]
) -> dict[str, float | tp.Callable]:
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


def estimate_half_life(
    durations: npt.NDArray[np.float64],
    observations: npt.NDArray[np.float64],
    replicate_indices: npt.NDArray[np.float64],
) -> dict[str, tp.Any]:
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


def estimate_half_life_per_replicate(
    durations: npt.NDArray[np.float64],
    observations: npt.NDArray[np.float64],
    replicate_indices: npt.NDArray[np.float64],
) -> dict[str, tp.Any]:
    metrics: dict[str, tp.Any] = {
        "decay_rate": [],
        "half_life": [],
        "half_life_std": [],
        "fit": [],
    }
    for repl in np.unique(replicate_indices):
        fit = estimate_half_life(
            durations=durations[replicate_indices == repl],
            observations=observations[replicate_indices == repl],
            replicate_indices=replicate_indices[replicate_indices == repl],
        )
        metrics["decay_rate"].append(fit["decay_rate"])
        metrics["half_life"].append(fit["half_life"])
        metrics["half_life_std"].append(fit["half_life_std"])
        metrics["fit"].append(fit["fit"])
    return metrics


def compute_cv_r2(
    durations: npt.NDArray[np.float64],
    observations: npt.NDArray[np.float64],
    replicate_indices: npt.NDArray[np.float64],
) -> float:
    r2s = []
    # estimate R2 goodness of fit score by running quasi-cross validation
    if len(np.unique(replicate_indices)) > 1:
        for repl_ind in np.unique(replicate_indices):
            mask = replicate_indices != repl_ind
            fit = fit_decay_curve(durations[mask], observations[mask])
            pred = tp.cast(tp.Callable, fit["fit"])(durations[~mask])
            r2s.append(r2_score(observations[~mask], pred))
    else:
        r2s.append(0)

    return np.array(r2s).mean()
