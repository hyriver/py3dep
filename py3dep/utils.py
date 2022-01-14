"""Utilities for Py3DEP."""
from typing import List, TypeVar

import numpy as np
import pygeoutils as geoutils
import rioxarray  # noqa: F401
import xarray as xr

try:
    import richdem as rd
except ImportError:
    rd = None

from .exceptions import MissingDependency

__all__ = ["deg2mpm", "fill_depressions"]
X = TypeVar("X", xr.DataArray, xr.Dataset)


def fill_depressions(dem_da: xr.DataArray) -> xr.DataArray:
    """Fill depressions and adjust flat areas in a DEM using `RichDEM <https://richdem.readthedocs.io>`__.

    Parameters
    ----------
    dem : xarray.DataArray or numpy.ndarray
        Digital Elevation Model.

    Returns
    -------
    xarray.DataArray
        Conditioned DEM after applying
        `depression filling <https://richdem.readthedocs.io/en/latest/depression_filling.html>`__
        and
        `flat area resolution <https://richdem.readthedocs.io/en/latest/flat_resolution.html>`__
        operations.
    """
    dem = dem_da.copy()
    if rd is None:
        raise MissingDependency

    nodata = -9999
    with xr.set_options(keep_attrs=True):  # type: ignore
        dem = dem.fillna(nodata)
        rda = rd.rdarray(dem, no_data=nodata)
        rda.projection = dem.rio.crs
        rda.geotransform = geoutils.transform2tuple(dem.rio.transform())
        rda = rd.FillDepressions(rda, epsilon=False)
        dem.data = rd.ResolveFlats(rda)
        return dem.where(dem > nodata, dem.attrs["nodatavals"][0])  # type: ignore


def deg2mpm(slope: xr.DataArray) -> xr.DataArray:
    """Convert slope from degrees to meter/meter.

    Parameters
    ----------
    slope : xarray.DataArray
        Slope in degrees.

    Returns
    -------
    xarray.DataArray
        Slope in meter/meter. The name is set to ``slope`` and the ``units`` attribute
        is set to ``m/m``.
    """
    with xr.set_options(keep_attrs=True):  # type: ignore
        nodata = slope._FillValue if hasattr(slope, "_FillValue") else slope.nodatavals[0]
        slope = slope.where(slope != nodata, drop=False)
        slope = np.tan(np.deg2rad(slope))
        slope.attrs["nodatavals"] = (np.nan,)
        if hasattr(slope, "_FillValue"):
            slope.attrs["_FillValue"] = np.nan
        slope.name = "slope"
        slope.attrs["units"] = "m/m"
    return slope


def rename_layers(ds: X, valid_layers: List[str]) -> X:
    """Rename layers in a dataset."""
    rename = {lyr: lyr.split(":")[-1].replace(" ", "_").lower() for lyr in valid_layers}
    rename.update({"3DEPElevation:None": "elevation"})

    if isinstance(ds, xr.DataArray):
        ds.name = rename[str(ds.name)]
    else:
        ds = ds.rename({n: rename[str(n)] for n in ds})

    return ds
