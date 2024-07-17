__version__ = "3.0.0"

import importlib.resources

with importlib.resources.path("illum", ".") as path:
    path = path.as_posix()

del importlib

from . import AngularPowerDistribution, MultiScaleData, SpectralPowerDistribution

MSD = MultiScaleData
APD = AngularPowerDistribution
SPD = SpectralPowerDistribution

from .UI import *
