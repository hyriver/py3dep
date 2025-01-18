"""Top-level package for Py3DEP."""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version

from py3dep import exceptions
from py3dep.geoops import deg2mpm, fill_depressions
from py3dep.print_versions import show_versions
from py3dep.py3dep import (
    add_elevation,
    check_3dep_availability,
    elevation_bycoords,
    elevation_bygrid,
    elevation_profile,
    get_dem,
    get_dem_vrt,
    get_map,
    query_3dep_sources,
    static_3dep_dem,
)

try:
    __version__ = version("py3dep")
except PackageNotFoundError:
    __version__ = "999"

__all__ = [
    "__version__",
    "add_elevation",
    "check_3dep_availability",
    "deg2mpm",
    "elevation_bycoords",
    "elevation_bygrid",
    "elevation_profile",
    "exceptions",
    "fill_depressions",
    "get_dem",
    "get_dem_vrt",
    "get_map",
    "query_3dep_sources",
    "show_versions",
    "static_3dep_dem",
]
