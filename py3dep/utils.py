"""Some utilities for Py3DEP."""
import os
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple, Union, ValuesView

import numpy as np
import rasterio as rio
import xarray as xr
from pygeoogc import MatchCRS
from rasterio import mask as rio_mask
from shapely.geometry import Polygon, box

from .exceptions import InvalidInputType


def create_dataset(
    content: bytes,
    geometry: Union[Polygon, Tuple[float, float, float, float]],
    geo_crs: str,
    name: str,
    fpath: Optional[Union[str, Path]] = None,
) -> Union[xr.Dataset, xr.DataArray]:
    """Create dataset from a response clipped by a geometry.

    Parameters
    ----------
    content : requests.Response
        The response to be processed
    geometry : Polygon or tuple of length 4
        The geometry for masking the data
    geo_crs : str
        The spatial reference of the input geometry
    name : str
        Variable name in the dataset
    fpath : str or Path, optinal
        The path save the file, defaults to None i.e., don't save as an image.

    Returns
    -------
    xarray.Dataset
        Generated xarray DataSet or DataArray
    """
    if not isinstance(geometry, (Polygon, tuple)):
        raise InvalidInputType("geometry", "Polygon or tuple of length 4")

    with rio.MemoryFile() as memfile:
        memfile.write(content)
        with memfile.open() as src:
            if src.nodata is None:
                try:
                    nodata = np.iinfo(src.dtypes[0]).max
                except ValueError:
                    nodata = np.nan
            else:
                nodata = np.dtype(src.dtypes[0]).type(src.nodata)

            if isinstance(geometry, Polygon):
                _geometry = [MatchCRS.geometry(geometry, geo_crs, src.crs)]
            else:
                _geometry = [box(*MatchCRS.bounds(geometry, geo_crs, src.crs))]

            masked, transform = rio_mask.mask(src, _geometry, crop=True, nodata=nodata)
            meta = src.meta
            meta.update(
                {
                    "width": masked.shape[2],
                    "height": masked.shape[1],
                    "transform": transform,
                    "nodata": nodata,
                }
            )

            if fpath is not None:
                check_dir(fpath)
                with rio.open(fpath, "w", **meta) as dest:
                    dest.write(masked)

            with rio.vrt.WarpedVRT(src, **meta) as vrt:
                ds = xr.open_rasterio(vrt)
                ds.data = masked
                try:
                    ds = ds.squeeze("band", drop=True)
                except ValueError:
                    pass
                ds.name = name

                ds.attrs["transform"] = transform
                ds.attrs["res"] = (transform[0], transform[4])
                ds.attrs["bounds"] = tuple(vrt.bounds)
                ds.attrs["nodatavals"] = vrt.nodatavals
                ds.attrs["crs"] = vrt.crs
    return ds


def wms_toxarray(
    r_dict: Dict[str, bytes],
    geometry: Union[Polygon, Tuple[float, float, float, float]],
    geo_crs: str,
    data_dir: Optional[Union[str, Path]] = None,
) -> Union[xr.DataArray, xr.Dataset]:
    """Convert responses from ``pygeoogc.wms_bybox`` to ``xarray.Dataset``.

    Parameters
    ----------
    r_dict : dict
        The output of ``wms_bybox`` function.
    geometry : Polygon or tuple
        The geometry to mask the data that should be in the same CRS as the r_dict.
    geo_crs : str
        The spatial reference of the input geometry.
    data_dir : str or Path, optional
        The directory to save the output as ``tiff`` images.

    Returns
    -------
    xarray.Dataset or xarray.DataAraay
        The dataset or data array based on the number of variables.
    """
    if data_dir is not None:
        check_dir(Path(data_dir, "dummy"))

    _fpath = {lyr: f'{lyr.split(":")[-1].lower().replace(" ", "_")}.tiff' for lyr in r_dict.keys()}
    fpath = {lyr: None if data_dir is None else Path(data_dir, fn) for lyr, fn in _fpath.items()}
    var_name = {lyr: f"{''.join(n for n in lyr.split('_')[:-1])}" for lyr in r_dict.keys()}

    ds = xr.merge(
        [
            create_dataset(r, geometry, geo_crs, var_name[lyr], fpath[lyr])
            for lyr, r in r_dict.items()
        ]
    )

    def mask_da(da):
        if not np.isnan(da.nodatavals[0]):
            msk = da < da.nodatavals[0] if da.nodatavals[0] > 0 else da > da.nodatavals[0]
            _da = da.where(msk, drop=True)
            _da.attrs["nodatavals"] = (np.nan,)
            return _da
        return da

    ds = ds.map(mask_da)

    if len(ds.variables) - len(ds.dims) == 1:
        ds = ds[list(ds.keys())[0]]
    return ds


def deg2mpm(da: xr.DataArray) -> xr.DataArray:
    """Convert ``xarray.Data[Array,set]`` from degree to meter/meter."""
    attrs = da.attrs
    da = np.tan(np.deg2rad(da))
    da.attrs = attrs
    da.name = "slope"
    da.attrs["units"] = "meters/meters"
    return da


def check_dir(
    fpath_itr: Optional[
        Union[ValuesView[Optional[Union[str, Path]]], List[Optional[Union[str, Path]]], str, Path]
    ]
) -> None:
    """Create parent directory for a file if doesn't exist."""
    if isinstance(fpath_itr, (str, Path)):
        fpath_itr = [fpath_itr]
    elif not isinstance(fpath_itr, Iterable):
        raise InvalidInputType("fpath_itr", "str or iterable")

    for f in fpath_itr:
        if f is None:
            continue

        parent = Path(f).parent
        if not parent.is_dir():
            try:
                os.makedirs(parent)
            except OSError:
                raise OSError(f"Parent directory cannot be created: {parent}")
