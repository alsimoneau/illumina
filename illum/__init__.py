__version__ = "2.2.2.20220620.23520770"

try:
    from . import MultiScaleData, pytools
    from .alternate import alternate
    from .batches import batches
    from .convert import convert
    from .domain import domain
    from .extract import extract
    from .failed import failed
    from .init import init
    from .inputs import inputs
    from .warp import warp
except ModuleNotFoundError:
    pass

try:
    from .compute import compute
except ModuleNotFoundError:
    print("Couldn't import Fortran functions")
