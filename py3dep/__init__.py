"""Top-level package for Py3DEP."""
from pkg_resources import DistributionNotFound, get_distribution

from .exceptions import InvalidInputType
from .print_versions import show_versions
from .py3dep import deg2mpm, elevation_bygrid, elevation_byloc, get_map

try:
    __version__ = get_distribution(__name__).version
except DistributionNotFound:
    # package is not installed
    __version__ = "999"
