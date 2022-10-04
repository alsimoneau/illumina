__version__ = "2.2.4.20221004.16174253"

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

import importlib.resources

with importlib.resources.path("illum", ".") as path:
    path = path.as_posix()

del importlib
