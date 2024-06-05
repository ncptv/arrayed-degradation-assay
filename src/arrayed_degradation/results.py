import io
import typing as tp
from collections import Counter
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image

from arrayed_degradation.decay_curve import fit_decay_curve
from arrayed_degradation.log import LOGGER
from arrayed_degradation.plots import (
    plot_summary_chart,
)

FLOAT_DATA_COLUMNS_ROUNDING = {
    "target_peak_areas": 2,
    "fraction_remaining": 3,
    "left_bound": 0,
    "right_bound": 0,
    "peak": 0,
    "prominence": 3,
}
FLOAT_RESULTS_COLUMNS_ROUNDING = {
    "decay_rate": 3,
    "decay_rate_std": 3,
    "half_life": 3,
    "half_life_std": 3,
    "r2_score": 3,
    "across_repl_corr": 4,
    "decay_rate_coef_of_var": 3,
    "half_life_coef_of_var": 3,
    "data_dissimilarity": 3,
}


TO_IGNORE_COLUMNS = ["epg_raw", "epg_clean", "epg_normed"]


def calculate_per_replicate_half_life(data: pd.DataFrame) -> pd.DataFrame:
    data = data.copy()
    data = data.loc[data["fraction_remaining"].notna()]
    per_replicate_results = data.groupby(["rna_id", "replicate"]).apply(
        lambda x: fit_decay_curve(x["timepoint"], x["fraction_remaining"])
    )
    return pd.DataFrame(
        per_replicate_results.tolist(), index=per_replicate_results.index
    ).reset_index()


def calculate_across_replicate_correlation(
    per_replicate_results: pd.DataFrame,
) -> float:
    """
    Calculate across replicate correlation. Correlation is calculated as follows:
    For each RNA and each replicate we have half life or decay rate estimated independently.
    We calculate correlation between all pairs of replicates, like this:
    (rna1_hf1, rna2_hf1, rna3_hf1, ...) vs (rna1_hf2, rna2_hf2m rna3_hf2, ...) etc.
    Then we average all these correlations.

    Args:
        per_replicate_results: DataFrame with columns: rna_id, replicate, half_life

    Returns:
        float: averaged corss-replicate rank correlation value
    """
    data = per_replicate_results.copy()
    all_one_repl_half_lifes = data.groupby("rna_id")["half_life"].apply(list).tolist()

    repl_count = map(len, all_one_repl_half_lifes)
    n_repl = Counter(repl_count).most_common(1)[0][
        0
    ]  # some replicates might have been dropped
    data = np.array([i for i in all_one_repl_half_lifes if len(i) == n_repl])
    corr = pd.DataFrame(data).corr("spearman").values
    corr = corr[np.tril_indices_from(corr, -1)]
    return np.nanmean(corr)


def calculate_coef_of_corr(
    per_replicate_results: pd.DataFrame, column: str
) -> dict[str, float]:
    return (
        per_replicate_results.groupby("rna_id")
        .apply(lambda x: x[column].std() / x[column].mean())
        .to_dict()
    )


def calculate_data_dissimilarity(data: pd.DataFrame) -> dict[str, float]:
    data = data.copy()
    data = data.loc[data["fraction_remaining"].notna() & data["failed"].isna()]
    scores = {}
    for rna_id, rna_data in data.groupby("rna_id"):
        dissimilarity = 0
        cnt = 0
        combs = combinations(rna_data["replicate"].unique(), 2)
        for rep1, rep2 in combs:
            rep1_data = rna_data.loc[rna_data["replicate"] == rep1]
            rep2_data = rna_data.loc[rna_data["replicate"] == rep2]
            common_timepoints = set(rep1_data.timepoint) & set(rep2_data.timepoint)
            cnt += len(common_timepoints)

            rep1_data = (
                rep1_data.loc[rep1_data.timepoint.isin(common_timepoints)][
                    "fraction_remaining"
                ]
                .astype(float)
                .values
            )
            rep2_data = (
                rep2_data.loc[rep2_data.timepoint.isin(common_timepoints)][
                    "fraction_remaining"
                ]
                .astype(float)
                .values
            )
            dissimilarity += np.abs((np.log(rep1_data) - np.log(rep2_data)) ** 2).sum()

        if cnt != 0:
            scores[rna_id] = dissimilarity / cnt
        else:
            scores[rna_id] = np.nan

    return scores


def process_plots(
    plots: list[Image.Image], all_results_df: pd.DataFrame, data_dir_path: Path
) -> None:
    results_file_name = "results.pdf"
    save_path = data_dir_path / results_file_name

    summary_plot = plot_summary_chart(all_results_df)
    plots.insert(0, summary_plot)
    LOGGER.info(f"Saving plots to {save_path}")

    with io.BytesIO() as output:
        plots[0].save(
            output,
            "PDF",
            resolution=200.0,
            save_all=True,
            append_images=plots[1:],
        )
        save_path.write_bytes(output.getvalue())


def process_results(
    results_list: dict[str, dict[str, tp.Any]], data_dir_path: Path
) -> None:
    all_results_df = (
        pd.DataFrame(results_list)
        .T.reset_index()
        .rename(columns={"index": "rna_id"})
        .astype({"rna_id": str})[
            [
                "rna_id",
                "decay_rate",
                "decay_rate_std",
                "half_life",
                "half_life_std",
                "r2_score",
            ]
        ]
    )
    all_data_df = pd.concat([i["data"] for i in results_list.values()])

    per_replicate_results = calculate_per_replicate_half_life(all_data_df)

    across_repl_corr = calculate_across_replicate_correlation(per_replicate_results)
    all_results_df["across_repl_corr"] = across_repl_corr
    decay_rate_coef_of_var = calculate_coef_of_corr(per_replicate_results, "decay_rate")
    all_results_df["decay_rate_coef_of_var"] = all_results_df.rna_id.map(
        decay_rate_coef_of_var
    )
    half_life_coef_of_var = calculate_coef_of_corr(per_replicate_results, "half_life")
    all_results_df["half_life_coef_of_var"] = all_results_df.rna_id.map(
        half_life_coef_of_var
    )
    data_dissimilarity = calculate_data_dissimilarity(all_data_df)
    all_results_df["data_dissimilarity"] = all_results_df.rna_id.map(data_dissimilarity)

    results_save_path = data_dir_path / "results.csv"
    LOGGER.info(f"Saving results to {results_save_path}")
    all_results_df[list(FLOAT_RESULTS_COLUMNS_ROUNDING.keys())] = (
        all_results_df[FLOAT_RESULTS_COLUMNS_ROUNDING.keys()]
        .astype(float)
        .round(FLOAT_RESULTS_COLUMNS_ROUNDING)
    )
    all_results_df.to_csv(results_save_path, index=False)

    data_save_path = data_dir_path / "data.csv"
    LOGGER.info(f"Saving data to {data_save_path}")
    all_data_df[list(FLOAT_DATA_COLUMNS_ROUNDING.keys())] = (
        all_data_df[FLOAT_DATA_COLUMNS_ROUNDING.keys()]
        .astype(float)
        .round(FLOAT_DATA_COLUMNS_ROUNDING)
    )
    all_data_df.drop(columns=TO_IGNORE_COLUMNS, errors="ignore").to_csv(
        data_save_path, index=False
    )

    plots = [r["plot"] for r in results_list.values()]
    process_plots(plots, all_results_df, data_dir_path)
