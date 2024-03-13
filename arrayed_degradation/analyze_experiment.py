import argparse
import typing as tp
from glob import glob
from pathlib import Path

from tqdm import tqdm

from analyze_single_rna import (
    analyze_single_rna,
)
from input_data import (
    process_input_data,
)
from log import LOGGER
from results import process_results


def analyze_experiment(
    data_dir_path: str,
    min_peak: float,
    max_peak: float,
    control_peak_min: float = 235,
    control_peak_max: float = 310,
    disable_control_peak: bool = False,
    remove_background: bool = False,
    width_param: float = 0.85,
) -> None:
    data_dir_path_obj = Path(data_dir_path)
    LOGGER.info(f"Analyzing experiment...")

    LOGGER.info("Reading input data...")
    all_data = process_input_data(data_dir_path_obj)

    LOGGER.info("Analyzing data...")
    all_results: dict[str, dict[str, tp.Any]] = {}
    for rna_id, rna_data in tqdm(all_data.groupby("rna_id")):
        try:
            rna_results = analyze_single_rna(
                rna_data,
                min_peak,
                max_peak,
                control_peak_min,
                control_peak_max,
                disable_control_peak,
                remove_background,
                width_param,
            )
        except Exception as e:
            LOGGER.error(f"Failed to analyze RNA {rna_id}. Error: {e}")
            continue
        all_results[rna_id] = rna_results

    LOGGER.info("Processing results...")
    process_results(all_results, data_dir_path_obj)


# pylint: disable=missing-docstring
def main(arguments: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        "--data_dir_path",
        help="""Path to the directory containing the data files. The directory should contain
        (1) one plate map file describing positioning of RNA samples and timepoints in hours for
        each replicate, (2) EPG files in FragmentAnalyzer format (one EPG file per replicate).
        Plate map file should be named 'plate_map.csv' and EPG files should be named in 'epg_0.csv',
        'epg_1.csv', etc. format.""",
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
        "--disable_control_peak",
        help=("Disable the p4p6 control peak used for normalization."),
        default=False,
        type=bool,
        action=argparse.BooleanOptionalAction,
    )
    parser.add_argument(
        "--remove_background",
        help="""Whether to remove background from target peak. If True it will draw a line at
        the base of detected peak and only the area above that line will be considered
        for half life calculations. It is not adviced to set it because of unstable
        behaviour and less repeatable results.""",
        default=False,
        type=bool,
        action=argparse.BooleanOptionalAction,
    )
    parser.add_argument(
        "--width_param",
        help="""Width param impacts significantly the width of the peak that is detected.
        The higher the value the wider bounds of the target peak. Do not exeed 0.95 as
        you will capture everything as peak. 0.75 is safe but quite conservative. 0.85
        captures mostly what you expect but can grab some noise.""",
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
        disable_control_peak=args.disable_control_peak,
        remove_background=args.remove_background,
        width_param=args.width_param,
    )


if __name__ == "__main__":
    main()
