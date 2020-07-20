"""Top-level package for Py3DEP."""

from .exceptions import InvalidInputType
from .py3dep import elevation_bygrid, elevation_byloc, get_map

__author__ = """Taher Chegini"""
__email__ = "cheginit@gmail.com"
__version__ = "0.1.3"

from ._version import get_versions

__version__ = get_versions()["version"]
del get_versions
