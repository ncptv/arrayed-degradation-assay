import typing as tp

import numpy as np
import pandas as pd

from decay_curve import (
    estimate_half_life,
    estimate_half_life_per_replicate,
)
from input_data import (
    EPG,
)
from peaks import (
    compute_peak_bounds,
    get_peak_area,
    normalize_trace_wrt_control,
)
from plots import (
    Plotter,
)


def analyze_single_rna(
    rna_data: pd.DataFrame,
    min_peak: float,
    max_peak: float,
    control_peak_min: float = 235,
    control_peak_max: float = 310,
    disable_control_peak: bool = False,
    remove_background: bool = False,
    width_param: float = 0.85,
) -> dict[str, tp.Any]:
    rna_data = rna_data.copy()
    plotter: Plotter = Plotter(
        sequence_name=rna_data.rna_id.values[0],
        timepoints=rna_data.timepoint.unique(),
        n_replicate=len(rna_data.replicate.unique()),
    )

    # normalize traces with respect to control peak so that the area under the control peak is equal to 1
    # for all timepoints and replicates
    if not disable_control_peak:
        rna_data["epg_normed"] = rna_data.epg_raw.apply(
            lambda x: normalize_trace_wrt_control(x, control_peak_min, control_peak_max)
        )

    # create arrays for storing data

    # timepoints; for example [0, 3, 6, 0, 3, 6, 0, 3, 6] means that there are three timepoints (0h, 3h, 6h) across 3 replicates
    timepoints = np.array([])

    # replicate indices for all timepoints; for example [0, 0, 0, 1, 1, 1, 2, 2, 2] means that the first
    # three timepoints are from replicate 0, the next three from replicate 1, and the last three from replicate 2
    replicate_indices = np.array([])

    # raw area under the target peak
    target_peak_areas = np.array([])

    # normalized area under the target peak with respect to t=0; for example [1, 0.5, 0.25, 1, 0.5, 0.25, 1, 0.5, 0.25]
    # means that for the first replicate at t_0 we had 100% of the sample, at t_1 we had 50%, and at t_2 we had 25%
    target_peak_areas_normed = np.array([])

    # iterate over replicates
    for replicate, replicate_data in rna_data.groupby("replicate"):
        # Compute bounds for the first peak (t=0); it will be used to compute peak areas for all remaning timepoints
        tp0_epg: EPG = replicate_data.loc[
            replicate_data.timepoint == 0, "epg_raw"
        ].values[0]
        left_bound, right_bound = compute_peak_bounds(
            tp0_epg.nucleotides,
            tp0_epg.trace,
            min_peak,
            max_peak,
            width_param,
        )

        this_replicate_peak_areas = []
        # iterate over timepoints of a given replicate
        for timepoint, timepoint_data in replicate_data.groupby("timepoint"):
            timepoints = np.append(timepoints, timepoint)
            replicate_indices = np.append(replicate_indices, replicate)

            epg: EPG = (
                timepoint_data.epg_normed.values[0]
                if not disable_control_peak
                else timepoint_data.epg_raw.values[0]
            )

            plotter.plot_sample_trace(
                epg.nucleotides,
                epg.trace,
                timepoint,
                left_bound,
                right_bound,
                control_peak_min if not disable_control_peak else None,
                control_peak_max if not disable_control_peak else None,
                replicate,
                remove_background,
            )

            # Compute area under the peak
            peak_area = get_peak_area(
                epg.nucleotides,
                epg.trace,
                left_bound,
                right_bound,
                remove_background,
            )
            this_replicate_peak_areas.append(peak_area)

        # example: [400, 200, 100] -> [1, 0.5, 0.25]
        this_replicate_normed_peak_areas = [
            peak_area / this_replicate_peak_areas[0]
            for peak_area in this_replicate_peak_areas
        ]

        target_peak_areas = np.append(target_peak_areas, this_replicate_peak_areas)
        target_peak_areas_normed = np.append(
            target_peak_areas_normed, this_replicate_normed_peak_areas
        )

    results = estimate_half_life(
        durations=timepoints,
        observations=target_peak_areas_normed,
        replicate_indices=replicate_indices,
    )

    plotter.combined_plotting(
        results,
        timepoints,
        np.array([]),
        target_peak_areas,
        target_peak_areas_normed,
    )

    results["plot"] = plotter.save_fig()

    return results
