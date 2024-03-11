import numpy as np
import numpy.typing as npt
from scipy.signal import find_peaks

from arrayed_degradation.exceptions import (
    PeakOutOfBoundsError,
)
from arrayed_degradation.input_data import EPG


def compute_peak_bounds(
    nucleotides: npt.NDArray[np.float64],
    trace: npt.NDArray[np.float64],
    min_peak: float,
    max_peak: float,
    width_param: float,
) -> tuple[float, float]:
    peaks, props = find_peaks(
        trace, prominence=(None, None), width=(None, None), rel_height=width_param
    )

    for p, _, l, r in sorted(
        zip(
            peaks,
            props["prominences"],
            props["left_ips"],
            props["right_ips"],
        ),
        key=lambda x: x[1],
        reverse=True,
    ):
        if nucleotides[p] >= min_peak and nucleotides[p] <= max_peak:
            left = l
            right = r
            break
    else:
        raise PeakOutOfBoundsError(
            f"Failed to find peak in provided bounds {min_peak}-{max_peak}"
        )

    return nucleotides[int(left)], nucleotides[int(right)]


def get_peak_area(
    nucleotides: npt.NDArray[np.float64],
    trace: npt.NDArray[np.float64],
    min_peak: float,
    max_peak: float,
    remove_background: bool,
) -> float:
    """
    Get peak area optionally normalized to the background.
    Area is calculated by integration of the intensity values over the nucleotide length.
    Integration is bound to the peak bounds.
    If remove_background is set to True the background area is subtracted from the peak area.
    Background area is defined as the area under the line connecting the first and last intensity value
    within the peak bounds. This is a simple way to estimate the background area and it is not perfect
    so it is advised not to remove background for further timepoints that might have migrating or split peaks.

    Args:
        nucleotides: npt.NDArray[np.float64]
            Array with sequence nucleotide length corresponding to each intensity reading.
        values: npt.NDArray[np.float64]
            Values of intensity of electropherogram.
        min_peak: float
            Left bound of peak. All values on the left will be masked and ignored.
        max_peak: float
            Right bound of peak. All values on the right will be masked and ignored.
        remove_background: bool
            Whether to remove background by drawing a line at the base of the peak
            (detemined by the left and right bound) and ignoring area that is underneath.
            It is advised not to remove background for target peak.

    Results:
        float:
            Nonnegative area of the peak.
    """
    peak_indices = np.where((nucleotides >= min_peak) & (nucleotides <= max_peak))

    nucleotides = nucleotides[peak_indices]
    values = trace[peak_indices]

    peak_area = np.trapz(y=values, x=nucleotides)

    if remove_background:
        background_area = np.trapz(
            y=[values[0], values[-1]], x=[nucleotides[0], nucleotides[-1]]
        )
        peak_area = peak_area - background_area
    return max(0, peak_area)


def normalize_trace_wrt_control(
    epg: EPG, control_peak_min: float, control_peak_max: float
) -> EPG:
    """Normalize trace of EPG with respect to the p4p6 control peak.
    After normalization the area under the control peak is 1, and all traces across replicates
    should be directly comparable.
    It is done to remove impact of pipetting errors and to make sure that the same amount of
    RNA is analyzed in each sample.
    """

    control_area = get_peak_area(
        epg.nucleotides,
        epg.trace,
        control_peak_min,
        control_peak_max,
        remove_background=True,
    )
    normed_trace = epg.trace / control_area
    return EPG(trace=normed_trace, nucleotides=epg.nucleotides)
