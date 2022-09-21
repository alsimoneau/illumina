__version__ = "2.2.4.20220921.15373229"

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

path = importlib.resources.files("illum").as_posix()
del importlib
