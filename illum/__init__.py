<<<<<<< HEAD
__version__ = "2.2.4.20221017.19280588"
=======
__version__ = "2.2.4.20221012.19375430"
>>>>>>> 2.2.4.20221012.19375430: Better multiple scale detection

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
