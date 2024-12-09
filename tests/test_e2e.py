from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd
from pytest import approx

from arrayed_degradation.analyze_experiment import (
    analyze_experiment,
)
from arrayed_degradation.tests.utils import (
    create_input_data,
)


def test_e2e():
    with TemporaryDirectory() as tmp_dir:
        create_input_data(tmp_dir)

        path = analyze_experiment(
            data_dir_path=tmp_dir,
            min_peak=1000,
            max_peak=2000,
            control_peak_min=100,
            control_peak_max=200,
            time_unit="h",
            disable_control_peak=False,
            remove_background=False,
            rel_height=0.75,
        )

        results = pd.read_csv(Path(path) / "results.csv", dtype={"rna_id": str})

    first_rna = results.loc[results.rna_id == "0001"]
    assert first_rna["half_life"].values[0] == approx(3, abs=1e-2)
    assert first_rna["decay_rate"].values[0] == approx(0.2308, abs=1e-2)
    assert first_rna["half_life_std"].values[0] == approx(0.0025, abs=1e-2)
    assert first_rna["r2_score"].values[0] == approx(0.999, abs=1e-2)

    second_rna = results.loc[results.rna_id == "0002"]
    assert second_rna["half_life"].values[0] == approx(1, abs=1e-2)
    assert second_rna["decay_rate"].values[0] == approx(0.691, abs=1e-2)
    assert second_rna["r2_score"].values[0] == approx(0.999, abs=1e-2)
