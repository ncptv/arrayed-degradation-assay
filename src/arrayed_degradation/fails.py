import pandas as pd


def fail_replicates_with_single_good_trace(rna_data: pd.DataFrame) -> pd.DataFrame:
    """Fail the whole replicate if there is only one good (not failed) trace in it.
    Reason for that is that having just one good trace in a replicate is pointless
    as we don't have anything to compare it to and we can't calculate the half life.

    Args:
        rna_data (pd.DataFrame): DataFrame with columns "replicate" and "failed".

    Returns:
        pd.DataFrame: DataFrame with failed column updated.
    """

    # if there is only one good trace in a replicate, we can fail the whole replicate
    replicates_with_single_good_trace_mask = (
        rna_data.groupby("replicate")["failed"].apply(lambda x: sum(x.isna())) == 1
    )
    replicates_to_fail = replicates_with_single_good_trace_mask.loc[
        replicates_with_single_good_trace_mask
    ].index
    rna_data.loc[
        rna_data.replicate.isin(replicates_to_fail) & rna_data["failed"].isna(),
        "failed",
    ] = "No other data points in this replicate."
    return rna_data


def fail_replicates_with_failed_tp0(rna_data: pd.DataFrame) -> pd.DataFrame:
    """Fail the whole replicate if timepoint 0 is failed.
    Reason for that is that we need the data at t=0 to calculate the half life.
    More specifically, we need the data at t=0 to calculate the initial amount
    of non degraded RNA and the fractions that are remaining at the other timepoints.

    Args:
        rna_data (pd.DataFrame): DataFrame with columns "replicate", "timepoint", "failed".

    Returns:
        pd.DataFrame:
    """

    # if there is no good data at t=0, we can fail the whole replicate
    tp0_failed_replicates = rna_data.loc[
        (rna_data.timepoint == 0) & (rna_data["failed"].notna()), "replicate"
    ].unique()
    rna_data.loc[
        (rna_data.replicate.isin(tp0_failed_replicates)) & (rna_data.failed.isna()),
        "failed",
    ] = "No data point at t=0."
    return rna_data


def fail_samples(
    rna_data: pd.DataFrame,
    message: str,
    replicate: int | list[int] | None = None,
    timepoint: float | int | list[float | int] | None = None,
) -> pd.DataFrame:
    """Fail sample or samples of given replicate (or replicates) and timepoint (or timepoints) with given message.
    Additionally fails the whole replicate if (1) there is only one good (not failed) trace in it,
    or (2) if timepoint 0 is failed.

    Args:
        rna_data (pd.DataFrame): DataFrame with columns "replicate", "timepoint", "failed".
        message (str): Message to put in the failed column.
        replicate (int | list[int] | None): Replicate or replicates to fail.
        timepoint (float | int | list[float  |  int] | None):  Timepoint or timepoints to fail.

    Returns:
        pd.DataFrame: DataFrame with failed column updated.
    """

    mask = rna_data.failed.isna()
    if replicate is not None:
        replicate = [replicate] if isinstance(replicate, int) else replicate
        mask &= rna_data.replicate.isin(replicate)
    if timepoint is not None:
        timepoint = (
            [timepoint]
            if isinstance(timepoint, float) or isinstance(timepoint, int)
            else timepoint
        )
        mask &= rna_data.timepoint.isin(timepoint)
    rna_data.loc[mask, "failed"] = message
    rna_data = fail_replicates_with_single_good_trace(rna_data)
    rna_data = fail_replicates_with_failed_tp0(rna_data)
    return rna_data
