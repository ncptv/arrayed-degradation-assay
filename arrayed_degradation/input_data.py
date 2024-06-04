import re
import string
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from arrayed_degradation.exceptions import EPGError
from arrayed_degradation.log import LOGGER

LABEL_FORMAT = r".+_TP\d+(?:\.\d+)?"


@dataclass
class EPG:
    """
    Class representing electropherogram data. Contains trace and nucleotides arrays (both of them must have equal size).

    Raises:
        EPGError: if trace and nucleotides arrays have different shapes
    """

    trace: np.ndarray
    nucleotides: np.ndarray

    def __post_init__(self) -> None:
        if self.trace.shape != self.nucleotides.shape:
            raise EPGError(
                f"Trace and nucleotides arrays should have the same shape. Got {self.trace.shape} and {self.nucleotides.shape}."
            )


def process_plate_map(plate_map: pd.DataFrame, time_unit: str) -> pd.DataFrame:
    """
    Process plate map from pivot format into tabular format.
    Decompose sample labels into df with RNA id, timepoint and well columns.
    Expected format of sample label is: `{rna_id}_TP{i}` where rna_id is
    a number determining id of RNA and i is a timepoint in hours.
    Input df should not contain any platemap specific column/row notation,
    these are assigned within that function.

    Example of processing.
    Input:
    0001_TP0  0001_TP0
    0001_TP3  0001_TP3
    Output:
    rna_id  timepoint  well   inplate_replicate
    0001    0          A1     0
    0001    0          A2     1
    0001    3          B1     0
    0001    3          B2     1
    """

    LOGGER.info("Processing plate map...")
    plate_map = plate_map.copy()

    # assign letters to rows and numbers to columns to obtain plate map specific positioning
    plate_map.index = [
        string.ascii_letters[i].upper() for i in range(plate_map.shape[0])
    ]
    plate_map.columns = range(1, plate_map.shape[1] + 1)
    plate_map.index.name = "row"  # type: ignore
    plate_map.columns.name = "column"  # type: ignore

    # unpivot plate map
    plate_map = (
        plate_map.stack(dropna=False).reset_index().rename(columns={0: "sample_id"})
    )
    plate_map["well"] = plate_map["row"] + plate_map["column"].astype(str)

    # checking label format
    plate_map["sample_id"] = plate_map["sample_id"].fillna("")
    proper_label = plate_map.sample_id.str.contains(LABEL_FORMAT, regex=True).fillna(
        False
    )

    if not proper_label.all():
        improper_labels = plate_map.loc[~proper_label, "well"].to_list()
        LOGGER.warning(
            f"""Empty/invalid sample labels detected on a plate map on these wells: {improper_labels}.
            Here are labels of these wells: {plate_map.loc[plate_map.well.isin(improper_labels), 'sample_id'].tolist()}.
            Check if the input data is correct and empty wells are intentional.
            Here is the correct format of the label: {{rna_id}}_TP{{i}} where rna_id is an id of RNA and i is a timepoint.
            Examples of correct label: `mRNA_0001_TP0`, `1234_TP3.5`."""
        )
        plate_map = plate_map.loc[proper_label]

    plate_map["timepoint"] = plate_map.sample_id.str.extract(
        r"TP(\d+(?:\.\d+)?)"
    ).astype(float)
    # convert timepoint to hours
    if time_unit == "m":
        plate_map["timepoint"] = plate_map["timepoint"] / 60

    plate_map["rna_id"] = plate_map.sample_id.str.extract(r"(.+)_TP")
    plate_map = plate_map.sort_values(["rna_id", "timepoint", "well"])
    plate_map.loc[:, "inplate_replicate"] = plate_map.groupby("sample_id").cumcount()

    return plate_map[["rna_id", "timepoint", "well", "inplate_replicate"]].reset_index(
        drop=True
    )


def process_epgs(epgs: pd.DataFrame) -> dict[str, EPG]:
    """
    Process Fragment Analyzer electropherogram data into EPG objects.
    First column in EPG df should contain size in nucleotides, next columns
    should contain electropherogram traces for each sample.
    Traces columns should contain {letter}{digit} well identifier like A6 or B12.
    """
    epgs = epgs.copy()
    well_ids = [re.search(r"([A-Z]+)(\d+)", i).group(0) for i in epgs.columns[1:]]  # type: ignore
    epgs.columns = ["Size (nt)"] + well_ids
    epgs = epgs.set_index("Size (nt)")
    return epgs.apply(
        lambda x: EPG(trace=x.values, nucleotides=x.index.values), axis=0
    ).to_dict()


def process_plate(plate_path: Path, time_unit: str) -> pd.DataFrame:
    """
    Process plate data from a single plate directory.

    Args:
        plate_path (Path):
            path to the plate directory
        time_unit (str):
            time unit used in the experiment; 'h' for hours, 'm' for minutes

    Raises:
        FileNotFoundError:
            raised if there is no `epg.csv` or `plate_map.csv` file in the plate directory

    Returns:
        pd.DataFrame:
            processed plate data along with EPGs assigned to each well
    """

    plate_map_path = plate_path / "plate_map.csv"
    epg_path = plate_path / "epg.csv"
    if not plate_map_path.exists():
        raise FileNotFoundError(f"Plate map file {plate_map_path} does not exist.")
    if not epg_path.exists():
        raise FileNotFoundError(f"EPG files {epg_path} do not exist.")

    # process plate map
    plate_map = pd.read_csv(plate_map_path, header=None)
    plate_map = process_plate_map(plate_map, time_unit)
    plate_map.loc[:, "plate_id"] = plate_path.name

    epgs = pd.read_csv(epg_path)
    epgs = process_epgs(epgs)

    return plate_map.assign(epg_raw=plate_map.apply(lambda x: epgs[x.well], axis=1))


def process_input_data(data_dir_path: Path, time_unit: str) -> pd.DataFrame:
    """Process input data from a directory containing plate subdirectories.

    Args:
        data_dir_path (Path):
            path to the directory containing plate subdirectories
        time_unit (str):
            time unit used in the experiment; 'h' for hours, 'm' for minutes

    Returns:
        pd.DataFrame:
            processed data from all plates in the data directory
    """

    if not data_dir_path.exists():
        raise FileNotFoundError(f"Data directory {data_dir_path} does not exist.")

    plate_paths = list(data_dir_path.glob("*"))
    plate_dir_paths = sorted([path for path in plate_paths if path.is_dir()])
    LOGGER.info(
        f"Found these plate subdirectories in the data directory: {[str(i) for i in plate_dir_paths]}"
    )

    plates_data = []
    for plate_path in plate_dir_paths:
        try:
            plates_data.append(process_plate(plate_path, time_unit))
        except FileNotFoundError as e:
            LOGGER.error(f"Error while processing plate {plate_path}: {e}")
            LOGGER.error(f"Skipping plate {plate_path}.")

    data = pd.concat(plates_data)
    data.loc[:, "replicate"] = data.groupby(["rna_id", "timepoint"])[
        ["plate_id", "inplate_replicate"]
    ].cumcount()
    data = (
        data[["rna_id", "timepoint", "replicate", "epg_raw"]]
        .sort_values(["rna_id", "timepoint", "replicate"])
        .reset_index(drop=True)
    )

    replicates = set(data.groupby("rna_id")["replicate"].nunique())
    n_rna = data.rna_id.nunique()
    timepoints = set(data.timepoint.unique())
    LOGGER.info(
        f"""Detected {n_rna} unique RNA ids, {replicates} replicates and {set(str(i) + "h" for i in timepoints)} timepoints."""
    )

    return data[["rna_id", "timepoint", "replicate", "epg_raw"]]
