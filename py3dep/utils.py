"""Utilities for Py3DEP."""
from __future__ import annotations

from typing import TYPE_CHECKING, Literal, Sequence, Union, overload

import numpy as np
import xarray as xr

if TYPE_CHECKING:
    import pyproj

    CRSTYPE = Union[int, str, pyproj.CRS]

from py3dep.exceptions import DependencyError

__all__ = ["deg2mpm", "fill_depressions"]


def fill_depressions(dem: xr.DataArray, outlets: Literal["min", "edge"] = "min") -> xr.DataArray:
    """Hydrocondition the DEM.

    This function uses `pyflwdir <https://deltares.github.io/pyflwdir/latest/>`__.
    It can be installed using ``pip install pyflwdir`` or
    ``conda install -c conda-forge pyflwdir``.

    Parameters
    ----------
    dem_da : xarray.DataArray or numpy.ndarray
        Digital Elevation Model.
    outlets : str, optional
        Initial basin outlet(s) at the edge of all cells
        (``edge``) or only the minimum elevation edge cell (``min``; default).

    Returns
    -------
    xarray.DataArray
        Conditioned DEM after applying ``fill_depressions`` function from
        `pyflwdir <https://deltares.github.io/pyflwdir/latest/>`__.
    """
    try:
        import pyflwdir
    except ImportError as ex:
        raise DependencyError from ex
    dem = dem.astype("f8")
    attrs = dem.attrs
    filled, _ = pyflwdir.dem.fill_depressions(dem.values, outlets=outlets, nodata=np.nan)
    dem = dem.copy(data=filled)
    dem.attrs = attrs
    dem = dem.rio.write_nodata(np.nan)
    return dem


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
    with xr.set_options(keep_attrs=True):
        if hasattr(slope, "_FillValue"):
            nodata = slope.attrs["_FillValue"]
        elif hasattr(slope, "nodatavals"):
            _nodata = slope.attrs["nodatavals"]
            nodata = _nodata[0] if isinstance(_nodata, Sequence) else _nodata
        else:
            nodata = np.nan
        slope = slope.where(slope != nodata, drop=False)

        def to_mpm(da: xr.DataArray) -> xr.DataArray:
            """Convert slope from degrees to meter/meter."""
            return xr.apply_ufunc(
                lambda x: np.tan(np.deg2rad(x)),
                da,
                vectorize=True,
            )

        slope = to_mpm(slope).compute()
        slope.attrs["nodatavals"] = (np.nan,)
        if hasattr(slope, "_FillValue"):
            slope.attrs["_FillValue"] = np.nan
        slope.name = "slope"
        slope.attrs["units"] = "m/m"
    return slope


@overload
def rename_layers(ds: xr.DataArray, valid_layers: list[str]) -> xr.DataArray:
    ...


@overload
def rename_layers(ds: xr.Dataset, valid_layers: list[str]) -> xr.Dataset:
    ...


def rename_layers(
    ds: xr.DataArray | xr.Dataset, valid_layers: list[str]
) -> xr.DataArray | xr.Dataset:
    """Rename layers in a dataset."""
    rename = {lyr: lyr.split(":")[-1].replace(" ", "_").lower() for lyr in valid_layers}
    rename.update({"3DEPElevation:None": "elevation"})

    if isinstance(ds, xr.DataArray):
        ds.name = rename[str(ds.name)]
    else:
        ds = ds.rename({n: rename[str(n)] for n in ds})
    return ds
