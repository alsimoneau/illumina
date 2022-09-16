__version__ = "2.2.4.20220916.15251599"

from . import AngularPowerDistribution, SpectralPowerDistribution
from .alternate import alternate
from .batches import batches
from .domain import domain
from .extract import extract
from .failed import failed
from .init import init
from .inputs import inputs
from .polar_warp import map_coordinates, polar_unwarp, polar_warp
from .warp import warp

try:
    from .compute import compute
except ModuleNotFoundError:
    print("Couldn't import Fortran functions")

import os as _os

try:
    path = _os.path.dirname(
        [
            path
            for path in _os.environ["PATH"].split(":")
            if path.endswith("illumina/bin")
        ][0]
    )
except IndexError:
    raise ValueError(
        "The 'illumina/bin' folder is not in the PATH environment variable."
    )

del _os
