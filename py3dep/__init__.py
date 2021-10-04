"""Top-level package for Py3DEP."""
from pkg_resources import DistributionNotFound, get_distribution

from .exceptions import (
    InvalidInputType,
    MissingAttribute,
    MissingColumns,
    MissingCRS,
    MissingDependency,
)
from .print_versions import show_versions
from .py3dep import elevation_bycoords, elevation_bygrid, get_map
from .utils import deg2mpm, fill_depressions

try:
    __version__ = get_distribution(__name__).version
except DistributionNotFound:
    __version__ = "999"

__all__ = [
    # Functions
    "fill_depressions",
    "get_map",
    "deg2mpm",
    "elevation_bycoords",
    "elevation_bygrid",
    "show_versions",
    # Exceptions
    "MissingColumns",
    "MissingCRS",
    "MissingDependency",
    "MissingAttribute",
    "InvalidInputType",
    # Constants
    "__version__",
]
