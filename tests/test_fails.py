import pandas as pd

from arrayed_degradation.fails import (
    fail_replicates_with_failed_tp0,
    fail_replicates_with_single_good_trace,
    fail_samples,
)


def test_fail_replicates_with_single_good_trace():
    rna_data = pd.DataFrame(
        {
            "replicate": [1, 1, 1, 2, 2, 2, 3, 3, 3],
            "timepoint": [0, 3, 6, 0, 3, 6, 0, 3, 6],
            "failed": [
                None,
                "some error",
                "some error",
                None,
                "some error",
                None,
                None,
                None,
                None,
            ],
        }
    )

    rna_data = fail_replicates_with_single_good_trace(rna_data)

    # Replicate 1 failed because it has only one good trace
    assert rna_data.loc[rna_data.replicate == 1, "failed"].isna().sum() == 0
    # Replicate 2 not since we have good 2 traces (and that's why the sum is 2)
    assert rna_data.loc[rna_data.replicate == 2, "failed"].isna().sum() == 2
    # And replicate 3 has no errors (so the sum is 3).
    assert rna_data.loc[rna_data.replicate == 3, "failed"].isna().sum() == 3


def testfail_replicates_with_failed_tp0():
    rna_data = pd.DataFrame(
        {
            "replicate": [1, 1, 1, 2, 2, 2, 3, 3, 3],
            "timepoint": [0, 3, 6, 0, 3, 6, 0, 3, 6],
            "failed": [
                "some error",
                None,
                None,
                None,
                "some error",
                None,
                "some error",
                None,
                "some error",
            ],
        }
    )

    rna_data = fail_replicates_with_failed_tp0(rna_data)

    # Replicate 1 failed because tp0 is failed so we have no reference and we fail entire replicate
    assert rna_data.loc[rna_data.replicate == 1, "failed"].isna().sum() == 0
    # Replicate 2 is not changed because tp=0 is good
    assert rna_data.loc[rna_data.replicate == 2, "failed"].isna().sum() == 2
    # Replicate 3 is failed because tp=0 is failed
    assert rna_data.loc[rna_data.replicate == 3, "failed"].isna().sum() == 0


def test_fail_samples():
    rna_data = pd.DataFrame(
        {
            "replicate": [1, 1, 1, 2, 2, 2, 3, 3, 3],
            "timepoint": [0, 3, 6, 0, 3, 6, 0, 3, 6],
            "failed": [
                None,
                None,
                None,
                None,
                "some error",
                None,
                None,
                None,
                None,
            ],
        }
    )

    # lets simulate failing some samples

    # we add single error that should not cause any additional errors
    rna_data = fail_samples(rna_data, "some error", replicate=1, timepoint=3)
    # 2 because we added 1 error
    assert rna_data.loc[rna_data.replicate == 1, "failed"].isna().sum() == 2
    # 2 because one sample was initially failed
    assert rna_data.loc[rna_data.replicate == 2, "failed"].isna().sum() == 2
    # 3 because all samples were initially good
    assert rna_data.loc[rna_data.replicate == 3, "failed"].isna().sum() == 3

    # we add error to all samples in replicate 2
    rna_data = fail_samples(rna_data, "some error", replicate=2)

    # repl 1 should not change
    assert rna_data.loc[rna_data.replicate == 1, "failed"].isna().sum() == 2
    # repl 2 should fail all samples
    assert rna_data.loc[rna_data.replicate == 2, "failed"].isna().sum() == 0
    # repl 3 should not change
    assert rna_data.loc[rna_data.replicate == 3, "failed"].isna().sum() == 3

    # we add error for timpoint 0 of last replicate that should cause the entire replicate to fail
    rna_data = fail_samples(rna_data, "some error", replicate=3, timepoint=0)

    # should not change
    assert rna_data.loc[rna_data.replicate == 1, "failed"].isna().sum() == 2
    # should not change
    assert rna_data.loc[rna_data.replicate == 2, "failed"].isna().sum() == 0
    # should fail all samples
    assert rna_data.loc[rna_data.replicate == 3, "failed"].isna().sum() == 0
