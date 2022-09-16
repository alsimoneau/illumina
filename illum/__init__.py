__version__ = "2.2.4.20220916.19375642"

import os

from . import (
    AngularPowerDistribution,
    MultiScaleData,
    PolarArray,
    SpectralPowerDistribution,
    utils,
)
from .UI import *

try:
    from .compute import compute
except ModuleNotFoundError:
    print("Couldn't import Fortran functions")


# Grab illumina code folder path from environement variables
try:
    path = os.path.dirname(
        [
            path
            for path in os.environ["PATH"].split(":")
            if path.endswith("illumina/bin")
        ][0]
    )
except IndexError:
    raise ValueError(
        "The 'illumina/bin' folder is not in the PATH environment variable."
    )

del os
