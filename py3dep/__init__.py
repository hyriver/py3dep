"""Top-level package for Py3DEP."""
from .exceptions import (
    InvalidInputType,
    InvalidInputValue,
    MissingColumns,
    MissingCRS,
    MissingDependency,
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
    import importlib.metadata
except ImportError:
    import importlib_metadata

    __version__ = importlib_metadata.version("py3dep")
else:
    __version__ = importlib.metadata.version("py3dep")

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
    "MissingColumns",
    "MissingCRS",
    "MissingDependency",
    "InvalidInputType",
    "InvalidInputValue",
    # Constants
    "__version__",
]
