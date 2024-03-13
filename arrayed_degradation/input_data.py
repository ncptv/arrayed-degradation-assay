import re
import string
from dataclasses import dataclass
from glob import glob
from pathlib import Path

import numpy as np
import pandas as pd

from log import LOGGER


@dataclass
class EPG:
    trace: np.ndarray
    nucleotides: np.ndarray


def process_plate_map(plate_map: pd.DataFrame) -> pd.DataFrame:
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
    plate_map["timepoint"] = plate_map.sample_id.str.extract(
        r"TP(\d+(?:\.\d+)?)"
    ).astype(float)
    plate_map["rna_id"] = plate_map.sample_id.str.extract(r"(.+)_TP")
    # drop empty wells
    empty_wells = plate_map.well.loc[
        plate_map[["timepoint", "rna_id"]].isna().any(axis=1)
    ].tolist()
    if empty_wells:
        LOGGER.warning(
            f"Empty/invalid sample labels detected on a plate map on these wells: {empty_wells}. "
            f"Here are labels of these wells: {plate_map.loc[plate_map.well.isin(empty_wells), 'sample_id'].tolist()} "
            "Check if the input data is correct and empty wells are intentional."
        )
        plate_map = plate_map.dropna(subset=["timepoint", "rna_id"])

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
    Traces columns should contain {letter}-{digit} well identifier like A6 or B12.
    """
    epgs = epgs.copy()
    well_ids = [re.search(r"([A-Z]+)(\d+)", i).group(0) for i in epgs.columns[1:]]  # type: ignore
    epgs.columns = ["Size (nt)"] + well_ids
    epgs = epgs.set_index("Size (nt)")
    return epgs.apply(
        lambda x: EPG(trace=x.values, nucleotides=x.index.values), axis=0
    ).to_dict()


def process_plate(plate_path: Path):
    plate_map_path = plate_path / "plate_map.csv"
    epg_path = plate_path / "epg.csv"
    if not plate_map_path.exists():
        raise FileNotFoundError(f"Plate map file {plate_map_path} does not exist.")
    if not epg_path.exists():
        raise FileNotFoundError(f"EPG files {epg_path} do not exist.")

    # process plate map
    plate_map = pd.read_csv(plate_map_path, header=None)
    plate_map = process_plate_map(plate_map)
    plate_map.loc[:, "plate_id"] = plate_path.name

    epgs = pd.read_csv(epg_path)
    epgs = process_epgs(epgs)

    return plate_map.assign(epg_raw=plate_map.apply(lambda x: epgs[x.well], axis=1))


def process_input_data(data_dir_path: Path) -> pd.DataFrame:
    if not data_dir_path.exists():
        raise FileNotFoundError(f"Data directory {data_dir_path} does not exist.")

    plate_paths = glob(str(data_dir_path / "*"))
    plate_dir_paths = [Path(i) for i in plate_paths if Path(i).is_dir()]
    LOGGER.info(f"Found these files in the data directory: {plate_dir_paths}")

    plates_data = []
    for plate_path in plate_dir_paths:
        plates_data.append(process_plate(Path(plate_path)))

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
