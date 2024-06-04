import argparse
import typing as tp

from pathy import Pathy
from tqdm import tqdm

from arrayed_degradation.analyze_single_rna import (
    analyze_single_rna,
)
from arrayed_degradation.exceptions import InVitroDegError
from arrayed_degradation.input_data import (
    process_input_data,
)
from arrayed_degradation.log import LOGGER
from arrayed_degradation.results import process_results


def analyze_experiment(
    data_dir_path: str,
    min_peak: float,
    max_peak: float,
    control_peak_min: float = 235,
    control_peak_max: float = 310,
    time_unit: str = "m",
    disable_control_peak: bool = False,
    remove_background: bool = False,
    rel_height: float = 0.85,
) -> str:
    data_dir_path_obj = Pathy.fluid(data_dir_path)
    LOGGER.info("Analyzing experiment...")

    LOGGER.info("Reading input data...")
    all_data = process_input_data(data_dir_path_obj, time_unit)

    LOGGER.info("Analyzing data...")
    all_results: dict[str, dict[str, tp.Any]] = {}

    sucesses, fails = 0, 0
    progress_bar = tqdm(total=all_data.rna_id.nunique())
    for i, (rna_id, rna_data) in enumerate(all_data.groupby("rna_id")):
        LOGGER.info(f"Analyzing RNA {rna_id}...")
        try:
            rna_results = analyze_single_rna(
                rna_data,
                min_peak,
                max_peak,
                control_peak_min,
                control_peak_max,
                disable_control_peak,
                remove_background,
                rel_height,
            )
            all_results[rna_id] = rna_results
            sucesses += 1
        except InVitroDegError as e:
            fails += 1
            LOGGER.error(f"Failed to analyze RNA {rna_id}. Error: {e}")
        if i > 0 and i % 10 == 0:
            progress_bar.set_postfix(
                {"rna_id": rna_id, "sucesses": sucesses, "fails": fails}, refresh=False
            )
            progress_bar.update(10)

    LOGGER.info(f"Analysis finished. Sucesses: {sucesses}, Fails: {fails}.")

    LOGGER.info("Processing results...")
    process_results(all_results, data_dir_path_obj)
    return data_dir_path


# pylint: disable=missing-docstring
def main(arguments: list[str] | None = None) -> str:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        "--data_dir_path",
        help="""Path to the directory containing the input data. The directory should contain
        subdirectories for each plate, each containing a plate_map.csv and epg.csv file.""",
        required=True,
        type=str,
    )
    parser.add_argument(
        "--min_peak",
        help="Minimum bound of initial peak, in nucleotides",
        required=True,
        type=float,
    )
    parser.add_argument(
        "--max_peak",
        help="Maximum bound of initial peak, in nucleotides",
        required=True,
        type=float,
    )
    parser.add_argument(
        "--control_peak_min",
        help="Minimum bound of control peak, in nucleotides",
        default=235,
        type=float,
    )
    parser.add_argument(
        "--control_peak_max",
        help="Maximum bound of control peak, in nucleotides",
        default=310,
        type=float,
    )
    parser.add_argument(
        "--time_unit",
        help="""Time unit for the timepoints in the plate_map.csv file.
        For hours use 'h', for minutes use 'm', for days 'd'.""",
        default="m",
        type=str,
    )
    parser.add_argument(
        "--disable_control_peak",
        help="Disable the control peak used for normalization.",
        default=False,
        type=bool,
        action=argparse.BooleanOptionalAction,
    )
    parser.add_argument(
        "--remove_background",
        help="""Whether to remove background from target peak. If True it will draw a line at
        the base of detected peak and only the area above that line will be considered
        for half life calculations. It is not adviced to set it because of unstable
        behaviour and less repeatable results for certain cases.""",
        default=False,
        type=bool,
        action=argparse.BooleanOptionalAction,
    )
    parser.add_argument(
        "--rel_height",
        help="""Relative height of the peak to be used for the `rel_height` bounding method.
        It determines the percentage of the peak height that will be used to draw a horizontal line
        to find the bounds of the peak. 
        Should be between 0 and 1 (although feasible values are between 0.5 and 0.95). The heigher the value the wider the bounds will be.""",
        default=0.85,
        type=float,
    )
    args = parser.parse_args(arguments)
    LOGGER.info(f"Analyzing experiment with the following args: {args}")
    return analyze_experiment(
        data_dir_path=args.data_dir_path,
        min_peak=args.min_peak,
        max_peak=args.max_peak,
        control_peak_min=args.control_peak_min,
        control_peak_max=args.control_peak_max,
        time_unit=args.time_unit,
        disable_control_peak=args.disable_control_peak,
        remove_background=args.remove_background,
        rel_height=args.rel_height,
    )


if __name__ == "__main__":
    main()
