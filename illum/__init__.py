__version__ = "2.2.4.20221021.16155395"

import importlib.resources

from . import AngularPowerDistribution as APD
from . import MultiScaleData as MSD
from . import PolarArray as PA
from . import SpectralPowerDistribution as SPD
from . import compute, utils
from .UI import *

with importlib.resources.path("illum", ".") as path:
    path = path.as_posix()

del importlib
