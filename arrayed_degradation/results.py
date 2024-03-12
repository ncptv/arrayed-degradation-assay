import typing as tp
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image

from log import LOGGER
from plots import (
    plot_summary_chart,
)


def calculate_across_replicate_correlation(
    all_one_repl_half_lifes: list[list[float]],
) -> float:
    """
    Calculate across replicate correlation. Correlation is calculated as follows:
    For each RNA and each replicate we have half life estimated independently.
    We calculate correlation between all pairs of replicates, like this:
    (rna1_hf1, rna2_hf1, rna3_hf1, ...) vs (rna1_hf2, rna2_hf2m rna3_hf2, ...) etc.
    Then we average all these correlations.

    Args:
        all_one_repl_half_lifes: list of lists of half lives for each RNA and each replicate.
            eg. [['rna1_hf1', 'rna1_hf2', 'rna1_hf3'], ['rna2_hf1', 'rna2_hf2', 'rna2_hf3'], ...]

    Returns:
        float: averaged corss-replicate rank correlation value
    """
    corr = pd.DataFrame(all_one_repl_half_lifes).corr("spearman").values
    corr = corr[np.tril_indices_from(corr, -1)]
    return np.nanmean(corr)


def process_plots(
    plots: list[Image.Image], all_results_df: pd.DataFrame, data_dir_path: Path
) -> None:
    save_path = data_dir_path / "results.pdf"
    summary_plot = plot_summary_chart(all_results_df)
    plots.insert(0, summary_plot)
    LOGGER.info(f"Saving plots to {save_path}")
    plots[0].save(
        save_path,
        "PDF",
        resolution=200.0,
        save_all=True,
        append_images=plots[1:],
    )


def process_results(
    results_list: dict[str, dict[str, tp.Any]], data_dir_path: Path
) -> None:
    all_one_repl_half_lives = [
        rna_results["one_repl_half_lifes"] for rna_results in results_list.values()
    ]
    across_replicate_correlation = calculate_across_replicate_correlation(
        all_one_repl_half_lives
    )

    all_results_df = (
        pd.DataFrame(results_list)
        .T.reset_index()
        .rename(columns={"index": "rna_id"})
        .assign(across_repl_corr=across_replicate_correlation)[
            [
                "rna_id",
                "decay_rate",
                "decay_rate_std",
                "half_life",
                "half_life_std",
                "r2_score",
                "across_repl_corr",
            ]
        ]
    )
    save_path = data_dir_path / "results.csv"
    LOGGER.info(f"Saving results to {save_path}")
    all_results_df.to_csv(save_path, index=False)

    plots = [r["plot"] for r in results_list.values()]
    process_plots(plots, all_results_df, data_dir_path)
