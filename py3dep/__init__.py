"""Top-level package for Py3DEP."""
from pkg_resources import DistributionNotFound, get_distribution

from .exceptions import InvalidInputType, InvalidInputValue, MissingColumns
from .print_versions import show_versions
from .py3dep import deg2mpm, elevation_bycoords, elevation_bygrid, get_map

try:
    __version__ = get_distribution(__name__).version
except DistributionNotFound:
    # package is not installed
    __version__ = "999"

__all__ = [
    # Functions
    "get_map",
    "deg2mpm",
    "elevation_bycoords",
    "elevation_bygrid",
    "show_versions",
    # Exceptions
    "InvalidInputType",
    "InvalidInputValue",
    "MissingColumns",
    # Constants
    "__version__",
]
