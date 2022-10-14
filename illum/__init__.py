__version__ = "2.2.4.20221014.15080532"

import importlib.resources

from . import (
    AngularPowerDistribution,
    MultiScaleData,
    PolarArray,
    SpectralPowerDistribution,
    compute,
    utils,
)
from .UI import *

with importlib.resources.path("illum", ".") as path:
    path = path.as_posix()

del importlib
