import re
import string
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from exceptions import LackOfTP0Error
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
    0001_TP0  0002_TP0
    0001_TP3  0002_TP3
    Output:
    rna_id  timepoint  well
    0001    0          A1
    0002    0          A2
    0001    3          B1
    0002    3          B2
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
    plate_map["rna_id"] = plate_map.sample_id.str.extract(r"(\d+)_TP")
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

    plate_map = plate_map.sort_values(["rna_id", "timepoint"])
    return plate_map[["rna_id", "timepoint", "well"]].reset_index(drop=True)


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


def process_input_data(plate_map_path: Path, epgs_paths: list[Path]) -> pd.DataFrame:
    if not plate_map_path.exists():
        raise FileNotFoundError(f"Plate map file not found: {plate_map_path}")
    for epgs_path in epgs_paths:
        if not epgs_path.exists():
            raise FileNotFoundError(f"EPG file not found: {epgs_path}")

    # process plate map
    plate_map = pd.read_csv(plate_map_path, header=None)
    plate_map = process_plate_map(plate_map)

    # process epgs; each EPG file should represent one replicate
    replicate_epgs: list[dict[str, EPG]] = []
    for epg_path in epgs_paths:
        epgs = pd.read_csv(epg_path)
        epgs = process_epgs(epgs)
        replicate_epgs.append(epgs)

    # merge plate map with epgs into one df
    # duplicate plate map for each replicate and concatenate
    data = pd.concat(
        [plate_map.assign(replicate=i) for i in range(len(replicate_epgs))]
    )
    # collect EPG for each well found on the plate map
    data["epg_raw"] = data.apply(lambda x: replicate_epgs[x.replicate][x.well], axis=1)
    data = data.sort_values(["rna_id", "replicate", "timepoint"]).reset_index(drop=True)

    n_replicate = data.replicate.nunique()
    n_rna = data.rna_id.nunique()
    timepoints = data.timepoint.unique()
    n_timepoints = len(timepoints)
    LOGGER.info(
        f"""Detected {n_rna} unique RNA ids, {n_replicate} replicates and {[str(i) + "h" for i in timepoints]} timepoints."""
    )
    expected_n_rows = n_rna * n_replicate * n_timepoints
    if expected_n_rows != data.shape[0]:
        LOGGER.warning(
            f"Detected {data.shape[0]} rows in the input data, but expected {n_rna} x {n_replicate} x {n_timepoints} "
            f"= {expected_n_rows} rows. Check if the input data is correct."
        )
    if data.loc[data.timepoint == 0].shape[0] != n_rna * n_replicate:
        raise LackOfTP0Error(
            f"Detected {data.loc[data.timepoint==0].shape[0]} rows in the input data for t=0, but expected {n_rna * n_replicate}. "
            "Check if the input data is correct; all RNAs should have timepoint 0 included. Aborting analysis."
        )

    return data[["rna_id", "replicate", "timepoint", "epg_raw"]]
