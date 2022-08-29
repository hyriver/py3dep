"""Top-level package for Py3DEP."""
from importlib.metadata import PackageNotFoundError, version

from .exceptions import (
    DependencyError,
    InputTypeError,
    InputValueError,
    MissingColumnError,
    MissingCRSError,
)
from .print_versions import show_versions
from .py3dep import (
    check_3dep_availability,
    elevation_bycoords,
    elevation_bygrid,
    elevation_profile,
    get_map,
    query_3dep_sources,
)
from .utils import deg2mpm, fill_depressions

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
