"""Top-level package for Py3DEP."""
from importlib.metadata import PackageNotFoundError, version

from packaging.version import Version

if Version(version("shapely")) > Version("1.9"):
    import os

    os.environ["USE_PYGEOS"] = "0"

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
    static_3dep_dem,
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
    "static_3dep_dem",
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
