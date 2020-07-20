"""Top-level package for Py3DEP."""
import pkg_resources

from .exceptions import InvalidInputType
from .py3dep import elevation_bygrid, elevation_byloc, get_map

try:
    __version__ = pkg_resources.get_distribution("xarray").version
except Exception:
    # Local copy or not installed with setuptools.
    # Disable minimum version checks on downstream libraries.
    __version__ = "999"
