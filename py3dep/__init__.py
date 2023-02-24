"""Top-level package for Py3DEP."""
from importlib.metadata import PackageNotFoundError, version

from py3dep.exceptions import (
    DependencyError,
    InputTypeError,
    InputValueError,
    MissingColumnError,
    MissingCRSError,
)
from py3dep.print_versions import show_versions
from py3dep.py3dep import (
    add_elevation,
    check_3dep_availability,
    elevation_bycoords,
    elevation_bygrid,
    elevation_profile,
    get_dem,
    get_map,
    query_3dep_sources,
    static_3dep_dem,
)
from py3dep.utils import deg2mpm, fill_depressions

try:
    __version__ = version("py3dep")
except PackageNotFoundError:
    __version__ = "999"

__all__ = [
    # Functions
    "fill_depressions",
    "get_map",
    "check_3dep_availability",
    "query_3dep_sources",
    "deg2mpm",
    "elevation_bycoords",
    "elevation_bygrid",
    "elevation_profile",
    "static_3dep_dem",
    "get_dem",
    "add_elevation",
    "show_versions",
    # Exceptions
    "MissingColumnError",
    "MissingCRSError",
    "DependencyError",
    "InputTypeError",
    "InputValueError",
    # Constants
    "__version__",
]
