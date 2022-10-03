__version__ = "2.2.4.20221003.19303390"

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
