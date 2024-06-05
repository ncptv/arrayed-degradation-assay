class InVitroDegError(Exception):
    pass


class PeakOutOfBoundsError(InVitroDegError):
    pass


class LackOfTP0Error(InVitroDegError):
    pass


class NormalizationError(InVitroDegError):
    pass


class EPGError(InVitroDegError):
    pass


class NoDataPointsError(InVitroDegError):
    pass
