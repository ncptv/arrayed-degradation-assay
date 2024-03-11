import typing as tp
from pathlib import Path

import pandas as pd
from PIL import Image

from log import LOGGER
from plots import plot_summary_chart


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
    all_results_df = (
        pd.DataFrame(results_list)
        .T.reset_index()
        .rename(columns={"index": "rna_id"})[
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
    save_path = data_dir_path / "results.csv"
    LOGGER.info(f"Saving results to {save_path}")
    all_results_df.to_csv(save_path, index=False)

    plots = [r["plot"] for r in results_list.values()]
    process_plots(plots, all_results_df, data_dir_path)
