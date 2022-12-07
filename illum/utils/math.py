# Mathematical functions
# Author: Alexandre Simoneau


import numpy as _np


def deg2mx(deg, lat):
    """Converts a longitudinal measurement from degrees to meters."""
    a = 6378137.0
    b = 6356752.3
    r = a * ((b / a * _np.tan(_np.deg2rad(lat))) ** 2 + 1) ** -0.5
    return deg * r


def deg2my(deg):
    """Converts a latitudinal measurement from degrees to meters."""
    return 40007863 * deg / 360


def safe_divide(a, b):
    """Safely divide two arrays, with 0 as a result of a division by 0."""
    with _np.errstate(divide="ignore", invalid="ignore"):
        c = _np.true_divide(a, b)
        c[c == _np.inf] = 0
        c = _np.nan_to_num(c)
    return c


def round_odd(n):
    return int(n - n % 2 + 1)
