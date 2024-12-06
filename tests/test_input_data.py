from itertools import product
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd
import pytest
from pytest import approx

from arrayed_degradation.exceptions import EPGError
from arrayed_degradation.input_data import (
    EPG,
    process_epgs,
    process_input_data,
    process_plate_map,
)
from arrayed_degradation.tests.utils import (
    create_input_data,
)


@pytest.fixture
def plate_map():
    return pd.DataFrame(
        [
            ["0001_TP0", "0001_TP3.5", "0001_TP6"],
            ["0001_TP0", "0001_TP3.5", "0001_TP6"],
            ["", "", ""],
            ["0002_TP0", "0002_TP3", "0002_TP6"],
            ["0002_TP0", "0002_TP3", "0002_TP6"],
        ]
    )


@pytest.fixture
def epgs():
    return pd.DataFrame(
        {
            "Size (nt)": [1, 2, 3, 4, 5],
        }
        | {
            f"{col}{row}": np.random.rand(5)
            for (col, row) in product("ABDE", range(1, 4))
        }
    )


def test_EPG():
    trace = np.array([1, 2, 3, 4, 5])
    nucleotides = np.array([0.1, 0.2, 0.3, 0.4, 0.5])

    epg = EPG(trace=trace, nucleotides=nucleotides)

    assert np.array_equal(epg.trace, trace)
    assert np.array_equal(epg.nucleotides, nucleotides)


def test_invalid_EPG():
    trace = np.array([1, 2, 3, 4, 5])
    nucleotides = np.array([0.1, 0.2, 0.3, 0.4])

    with pytest.raises(EPGError):
        EPG(trace=trace, nucleotides=nucleotides)


def test_process_plate_map(plate_map):
    processed_plate_map = process_plate_map(plate_map, "h")

    assert processed_plate_map.shape == (12, 4)
    assert processed_plate_map["rna_id"].unique().tolist() == ["0001", "0002"]
    assert sorted(processed_plate_map["timepoint"].unique().tolist()) == approx(
        [0, 3, 3.5, 6]
    )
    picked_well = processed_plate_map.loc[processed_plate_map.well == "E1", :]
    assert picked_well["rna_id"].values[0] == "0002"
    assert picked_well["timepoint"].values[0] == 0
    assert picked_well["inplate_replicate"].values[0] == 1


def test_process_epgs(epgs):
    processed_epg = process_epgs(epgs)
    assert len(processed_epg) == 12
    assert isinstance(processed_epg["A1"], EPG)
    assert all([len(i.trace) == 5 for i in processed_epg.values()])
    assert all([len(i.nucleotides) == 5 for i in processed_epg.values()])
    assert all(
        [(i.nucleotides == [1, 2, 3, 4, 5]).all() for i in processed_epg.values()]
    )


def test_process_input_data():
    with TemporaryDirectory() as tmpdir:
        create_input_data(
            tmpdir, seqs=["0001", "0002"], tps=[0, 3, 6], n_repl=2, half_lives=[3, 1]
        )
        data = process_input_data(Path(tmpdir), "h")
    assert data.shape == (12, 4)
    assert data["rna_id"].unique().tolist() == ["0001", "0002"]
    assert sorted(data["timepoint"].unique().tolist()) == approx([0, 3, 6])
