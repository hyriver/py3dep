"""Utilities for Py3DEP."""
import tempfile
import uuid
from pathlib import Path
from typing import Dict, List, TypeVar

import numpy as np
import pygeoutils as geoutils
import rasterio as rio
import rasterio.warp as rio_warp
import rioxarray as rxr
import xarray as xr

try:
    import richdem as rd
except ImportError:
    rd = None

from .exceptions import MissingAttribute, MissingDependency

__all__ = ["deg2mpm", "fill_depressions"]
X = TypeVar("X", xr.DataArray, xr.Dataset)


def reproject_gtiff(
    r_dict: Dict[str, bytes],
    crs: str,
) -> xr.DataArray:
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
    xarray.DataArray
        Reprojected data array.
    """
    tmp_dir = tempfile.gettempdir()
    var_name = {lyr: "_".join(lyr.split("_")[:-3]) for lyr in r_dict.keys()}
    attrs = geoutils.get_gtiff_attrs(next(iter(r_dict.values())), ("y", "x"))

    def to_dataset(lyr: str, resp: bytes) -> Path:
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
                    ds = rxr.open_rasterio(vrt)
                    ds = ds.squeeze("band", drop=True)
                    ds = ds.sortby(attrs.dims[0], ascending=False)
                    ds.attrs["crs"] = attrs.crs.to_string()
                    ds.attrs["transform"] = attrs.transform
                    ds.attrs["nodatavals"] = (attrs.nodata,)
                    ds.name = var_name[lyr]
                    fpath = Path(tmp_dir, f"{uuid.uuid4().hex}.nc")
                    ds.to_netcdf(fpath)
                    return fpath

    ds: xr.Dataset = xr.open_mfdataset(  # type: ignore
        (to_dataset(lyr, resp) for lyr, resp in r_dict.items()),
        parallel=True,
        decode_coords="all",
    )
    da: xr.DataArray = ds[list(ds.keys())[0]]
    da.attrs["crs"] = crs
    da.attrs["nodatavals"] = (attrs.nodata,)
    _transform = geoutils.get_transform(da, attrs.dims)[0]
    transform = tuple(getattr(_transform, c) for c in ["a", "b", "c", "d", "e", "f"])
    da = da.sortby(attrs.dims[0], ascending=False)
    da.attrs["transform"] = transform
    da.attrs["res"] = (_transform.a, _transform.e)

    for attr in ("scales", "offsets"):
        if attr in da.attrs and not isinstance(da.attrs[attr], tuple):
            da.attrs[attr] = (da.attrs[attr],)
    return da


def fill_depressions(
    dem: xr.DataArray,
) -> xr.DataArray:
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
    attrs = slope.attrs
    slope = np.tan(np.deg2rad(slope))
    slope.attrs = attrs
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
        ds = ds.rename({n: rename[n] for n in ds.keys()})

    return ds
