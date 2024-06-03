import typing as tp

import numpy as np
import pandas as pd

from arrayed_degradation.decay_curve import (
    estimate_half_life,
)
from arrayed_degradation.exceptions import (
    NoDataPointsError,
    NormalizationError,
    PeakOutOfBoundsError,
)
from arrayed_degradation.fails import fail_samples
from arrayed_degradation.input_data import (
    EPG,
)
from arrayed_degradation.log import LOGGER
from arrayed_degradation.peaks import (
    compute_peak_bounds,
    get_peak_area,
    normalize_trace_wrt_control,
)
from arrayed_degradation.plots import (
    Plotter,
)


def normalize_traces(
    rna_data: pd.DataFrame,
    control_peak_min: float,
    control_peak_max: float,
) -> pd.DataFrame:
    rna_data["epg_normed"] = None
    for repl, data_repl in rna_data.groupby("replicate"):
        for timepoint, data_tp in data_repl.groupby("timepoint"):
            try:
                normed_epg = normalize_trace_wrt_control(
                    data_tp.epg_raw.values[0], control_peak_min, control_peak_max
                )
                rna_data.loc[
                    (rna_data.replicate == repl) & (rna_data.timepoint == timepoint),
                    "epg_normed",
                ] = normed_epg
            except NormalizationError as e:
                rna_data = fail_samples(
                    rna_data,
                    "Normalization error",
                    replicate=repl,
                    timepoint=timepoint,
                )
                LOGGER.warn(
                    f"Failed to normalize trace for replicate {repl + 1}, timepoint {timepoint}. Failing datapoint. Error: {e}"
                )
                if timepoint == 0:
                    LOGGER.warn(
                        "Dropping remaming timepoints for this replicate as there is no correct tp=0."
                    )
                    break

    return rna_data


def analyze_single_rna(
    rna_data: pd.DataFrame,
    min_peak: float,
    max_peak: float,
    control_peak_min: float = 235,
    control_peak_max: float = 310,
    disable_control_peak: bool = False,
    remove_background: bool = False,
    rel_height: float = 0.85,
) -> dict[str, tp.Any]:
    rna_data = rna_data.copy()
    rna_data["failed"] = None

    # normalize traces with respect to control peak so that the area under the control peak is equal to 1
    # for all timepoints and replicates
    if not disable_control_peak:
        rna_data = normalize_traces(rna_data, control_peak_min, control_peak_max)
        if rna_data.failed.notna().all():
            raise NoDataPointsError(
                "All data points were dropped due to normalization errors. Check the logs for more information."
            )

    # raw area under the target peak
    rna_data["target_peak_areas"] = None

    # normalized area under the target peak with respect to t=0; for example [1, 0.5, 0.25, 1, 0.5, 0.25, 1, 0.5, 0.25]
    # means that for the first replicate at t_0 we had 100% of the sample, at t_1 we had 50%, and at t_2 we had 25%
    rna_data["fraction_remaining"] = None

    rna_data[["left_bound", "right_bound", "peak", "prominence"]] = None

    # iterate over replicates
    failed_mask = rna_data.failed.notna()
    for replicate, replicate_data in rna_data.loc[~failed_mask].groupby("replicate"):
        replicate_mask = (rna_data.replicate == replicate) & (rna_data.failed.isna())
        tp0_left_bound: float | None = None
        tp0_right_bound: float | None = None
        # iterate over timepoints of a given replicate
        for timepoint, timepoint_data in replicate_data.groupby("timepoint"):
            timepoint_mask = replicate_mask & (rna_data.timepoint == timepoint)
            sample_mask = replicate_mask & timepoint_mask
            epg: EPG = (
                timepoint_data.epg_normed.values[0]
                if not disable_control_peak
                else timepoint_data.epg_raw.values[0]
            )
            try:
                left_bound, peak, right_bound, prominence = compute_peak_bounds(
                    epg.nucleotides,
                    epg.trace,
                    min_peak if tp0_left_bound is None else tp0_left_bound - 200,
                    max_peak if tp0_right_bound is None else tp0_right_bound + 200,
                    rel_height=rel_height,
                )
            except PeakOutOfBoundsError as e:
                rna_data = fail_samples(
                    rna_data,
                    "Peak bounds error",
                    replicate=replicate,
                    timepoint=timepoint,
                )
                LOGGER.warn(
                    f"Failed to compute peak bounds for replicate {replicate + 1}. Error: {e}"
                )
                continue

            rna_data.loc[sample_mask, "left_bound"] = left_bound
            rna_data.loc[sample_mask, "right_bound"] = right_bound
            rna_data.loc[sample_mask, "peak"] = peak
            rna_data.loc[sample_mask, "prominence"] = prominence
            if tp0_left_bound is None or tp0_right_bound is None:
                tp0_left_bound = left_bound
                tp0_right_bound = right_bound

            if tp0_left_bound > peak or peak > tp0_right_bound:
                LOGGER.warn(
                    f"Peak for replicate {replicate + 1}, timepoint {timepoint} are not within the bounds of the peak at t=0."
                )
                rna_data = fail_samples(
                    rna_data,
                    "Peak shift error",
                    replicate=replicate,
                    timepoint=timepoint,
                )
                continue

            # Compute area under the peak
            peak_area = get_peak_area(
                epg.nucleotides,
                epg.trace,
                tp0_left_bound,
                tp0_right_bound,
                remove_background,
                prominence,
            )
            rna_data.loc[sample_mask, "target_peak_areas"] = peak_area
            if peak_area < 1e-2:
                LOGGER.warn(
                    f"Peak area for replicate {replicate + 1}, timepoint {timepoint} is 0."
                )
                rna_data = fail_samples(
                    rna_data,
                    "Peak area error",
                    replicate=replicate,
                    timepoint=timepoint,
                )

        else:
            # example: [400, 200, 100] -> [1, 0.5, 0.25]
            replicate_mask = (rna_data.replicate == replicate) & (
                rna_data.failed.isna()
            )
            this_replicate_peak_areas = rna_data.loc[
                replicate_mask, "target_peak_areas"
            ].values
            this_replicate_normed_peak_areas = [
                peak_area / this_replicate_peak_areas[0]
                for peak_area in this_replicate_peak_areas
            ]
            rna_data.loc[
                replicate_mask,
                "fraction_remaining",
            ] = this_replicate_normed_peak_areas

            if (np.diff(this_replicate_normed_peak_areas) > 0).any():
                LOGGER.warn(
                    f"Peak areas for replicate {replicate + 1} are not monotonically decreasing."
                )
                rna_data = fail_samples(
                    rna_data,
                    "Peak areas are not decreasing",
                    replicate=replicate,
                )

    failed_mask = rna_data.failed.notna()
    if failed_mask.all():
        raise NoDataPointsError(
            "All data points were dropped due to errors during processing."
        )

    results = estimate_half_life(
        durations=rna_data.loc[~failed_mask, "timepoint"].values,
        observations=rna_data.loc[~failed_mask, "fraction_remaining"].values,
        replicate_indices=rna_data.loc[~failed_mask, "replicate"].values,
    )

    plotter: Plotter = Plotter(
        rna_data=rna_data,
        fit_results=results,
        remove_background=remove_background,
        control_peak_min=control_peak_min if not disable_control_peak else None,
        control_peak_max=control_peak_max if not disable_control_peak else None,
    )
    plotter.combined_plotting()
    results["plot"] = plotter.save_fig()
    results["data"] = rna_data

    return results
