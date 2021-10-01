"""Utilities for Py3DEP."""
import tempfile
import uuid
from pathlib import Path
from typing import Dict, Union

import numpy as np
import pygeoutils as geoutils
import rasterio as rio
import rasterio.warp as rio_warp
import xarray as xr

try:
    import richdem as rd
except ImportError:
    rd = None

from .exceptions import MissingAttribute, MissingDependency

__all__ = ["deg2mpm", "fill_depressions"]


def reproject_gtiff(
    r_dict: Dict[str, bytes],
    crs: str,
) -> Union[xr.DataArray, xr.Dataset]:
    """Reproject a GTiff response into another CRS.

    Parameters
    ----------
    r_dict : dict
        Dictionary of (Geo)Tiff byte responses where keys are some names that are used
        for naming each responses, and values are bytes.
    crs : str
        The target spatial reference system.

    Returns
    -------
    xarray.DataArray or xarray.Dataset
        Reprojected data array.
    """
    tmp_dir = tempfile.gettempdir()
    var_name = {lyr: "_".join(lyr.split("_")[:-3]) for lyr in r_dict.keys()}
    attrs = geoutils.get_gtiff_attrs(next(iter(r_dict.values())), ("y", "x"))

    def to_dataset(lyr: str, resp: bytes) -> xr.DataArray:
        with rio.MemoryFile() as memfile:
            memfile.write(resp)
            with memfile.open() as src:
                transform, width, height = rio_warp.calculate_default_transform(
                    src.crs, crs, src.width, src.height, *src.bounds
                )
                kwargs = src.meta.copy()
                kwargs.update(
                    {"crs": crs, "transform": transform, "width": width, "height": height}
                )

                with rio.vrt.WarpedVRT(src, **kwargs) as vrt:
                    if crs != src.crs:
                        for i in range(1, src.count + 1):
                            rio_warp.reproject(
                                source=rio.band(src, i),
                                destination=rio.band(vrt, i),
                                src_transform=src.transform,
                                src_crs=src.crs,
                                dst_transform=transform,
                                crs=crs,
                                resampling=rio_warp.Resampling.bilinear,
                            )
                    ds = xr.open_rasterio(vrt)
                    ds = ds.squeeze("band", drop=True)
                    ds = ds.sortby(attrs.dims[0], ascending=False)
                    ds.attrs["crs"] = attrs.crs.to_string()
                    ds.attrs["transform"] = attrs.transform
                    ds.name = var_name[lyr]
                    fpath = Path(tmp_dir, f"{uuid.uuid4().hex}.nc")
                    ds.to_netcdf(fpath)
                    return fpath

    ds = xr.open_mfdataset(
        (to_dataset(lyr, resp) for lyr, resp in r_dict.items()),
        parallel=True,
    )
    ds = ds[list(ds.keys())[0]]
    ds.name = "elevation"
    ds.attrs["units"] = "meters"
    ds.attrs["crs"] = crs
    ds.attrs["nodatavals"] = (attrs.nodata,)
    _transform = geoutils.get_transform(ds, attrs.dims)[0]
    transform = tuple(getattr(_transform, c) for c in ["a", "b", "c", "d", "e", "f"])
    ds = ds.sortby(attrs.dims[0], ascending=False)
    ds.attrs["transform"] = transform
    ds.attrs["res"] = (_transform.a, _transform.e)

    for attr in ("scales", "offsets"):
        if attr in ds.attrs and not isinstance(ds.attrs[attr], tuple):
            ds.attrs[attr] = (ds.attrs[attr],)
    return ds


def fill_depressions(
    dem: Union[xr.DataArray, xr.Dataset],
) -> xr.DataArray:
    """Fill depressions and adjust flat areas in a DEM using RichDEM.

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
    if rd is None:
        raise MissingDependency

    if not hasattr(dem, "nodatavals"):
        raise MissingAttribute("dem", "nodatavals")

    no_data = dem.nodatavals[0] if isinstance(dem.nodatavals, tuple) else dem.nodatavals
    rda = rd.rdarray(dem, no_data=no_data)
    rda.projection = dem.crs
    rda.geotransform = dem.transform
    rda = rd.FillDepressions(rda, epsilon=False)
    dem.data = rd.ResolveFlats(rda)
    return dem


def deg2mpm(da: xr.DataArray) -> xr.DataArray:
    """Convert slope from degree to meter/meter.

    Parameters
    ----------
    da : xarray.DataArray
        Slope in degree.

    Returns
    -------
    xarray.DataArray
        Slope in meter/meter. The name is set to ``slope`` and the ``units`` attribute
        is set to ``m/m``.
    """
    attrs = da.attrs
    da = np.tan(np.deg2rad(da))
    da.attrs = attrs
    da.name = "slope"
    da.attrs["units"] = "m/m"
    return da
