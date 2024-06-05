import numpy as np
import numpy.typing as npt
from scipy.signal import find_peaks

from arrayed_degradation.exceptions import (
    NormalizationError,
    PeakOutOfBoundsError,
)
from arrayed_degradation.input_data import EPG


def compute_peak_bounds(
    nucleotides: npt.NDArray[np.float64],
    trace: npt.NDArray[np.float64],
    min_peak: float,
    max_peak: float,
    rel_height: float,
) -> tuple[float, float, float, float]:
    """Compute the bounds of the peak within the provided nucleotide range.

    This function takes in the nucleotide values and the corresponding trace values,
    and computes the bounds of the peak within the specified range. The peak bounds
    are determined based on the width parameter.

    Args:
        nucleotides (npt.NDArray[np.float64]):
            An array of nucleotide values.
        trace (npt.NDArray[np.float64]):
            An array of trace values.
        min_peak (float):
            The minimum value of the peak range.
        max_peak (float):
            The maximum value of the peak range.
        rel_height (float):
            The relative height parameter used to determine the peak width.

    Raises:
        PeakOutOfBoundsError: Raised when a peak within the specified bounds is not found.

    Returns:
        tuple[float, float, float, float]: A tuple containing the left bound, peak, right bound, and prominance of the peak.
    """

    peaks, props = find_peaks(
        trace, prominence=(None, None), width=(None, None), rel_height=rel_height
    )

    for p, pr, l, r in sorted(
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
            peak = p
            prominance = pr
            break
    else:
        raise PeakOutOfBoundsError(
            f"Failed to find peak in provided bounds {min_peak}-{max_peak}"
        )

    return (
        nucleotides[int(left)],
        nucleotides[int(peak)],
        nucleotides[int(right)],
        prominance,
    )


def get_peak_area(
    nucleotides: npt.NDArray[np.float64],
    trace: npt.NDArray[np.float64],
    min_peak: float,
    max_peak: float,
    remove_background: bool,
    prominence: float | None = None,
) -> float:
    """
    Get peak area optionally normalized to the background.
    Area is calculated by integration of the intensity values over the nucleotide length.
    Integration is bound to the peak bounds.
    There are two ways of background removal. The first one is to draw a line at the base of the peak
    (detemined by the left and right bound) and ignore area that is underneath. The area is trapezoidal.
    The second one is to draw a line at the height of the peak minus the prominence and ignore area that is
    underneath. The area is rectangular. The second method is more robust to the peak shape and noise.

    Args:
        nucleotides: npt.NDArray[np.float64]
            Array with sequence nucleotide length corresponding to each intensity reading.
        trace: npt.NDArray[np.float64]
            Values of intensity of electropherogram.
        min_peak: float
            Left bound of peak. All values on the left will be masked and ignored.
        max_peak: float
            Right bound of peak. All values on the right will be masked and ignored.
        remove_background: bool
            Whether to remove background by drawing a line at the base of the peak
            (detemined by the left and right bound) and ignoring area that is underneath.
            It is advised not to remove background for target peak.
        prominence: float | None
            Prominence of the peak. If not None and remove_background is set to True the background area
            is calculated as the area under the line drawn at the height of the peak minus the prominence
            (so it is rectangular). If None the background area is calculated as the area under the line
            connecting the first and last intensity value within the peak bounds.

    Returns:
        float:
            Nonnegative area of the peak.
    """
    peak_indices = np.where((nucleotides >= min_peak) & (nucleotides <= max_peak))

    nucleotides = nucleotides[peak_indices]
    values = trace[peak_indices]

    peak_area = np.trapz(y=values, x=nucleotides)

    if remove_background:
        if prominence is None:
            background_area = np.trapz(
                y=[values[0], values[-1]], x=[nucleotides[0], nucleotides[-1]]
            )
        else:
            max_value = max(values)
            base_level = max_value - prominence
            background_area = (max_peak - min_peak) * base_level
        peak_area = peak_area - background_area
    return max(0, peak_area)


def normalize_trace_wrt_control(
    epg: EPG, control_peak_min: float, control_peak_max: float
) -> EPG:
    """Normalize trace of EPG with respect to the control peak.
    After normalization the area under the control peak is 1, and all traces across replicates
    should be directly comparable.
    It is done to remove impact of pipetting errors and to make sure that the same amount of
    RNA is analyzed in each sample.

    Args:
        epg: EPG
            Electropherogram to be normalized.
        control_peak_min: float
            Left bound of control peak.
        control_peak_max: float
            Right bound of control peak.

    Returns:
        EPG:
            New instance of normalized electropherogram.
    """

    control_area = get_peak_area(
        epg.nucleotides,
        epg.trace,
        control_peak_min,
        control_peak_max,
        remove_background=True,
    )
    if control_area == 0:
        raise NormalizationError(
            "Control peak area is 0. Check the bounds of the control peak."
        )
    normed_trace = epg.trace / control_area
    return EPG(trace=normed_trace, nucleotides=epg.nucleotides)
