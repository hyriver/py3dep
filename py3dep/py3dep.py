"""Get data from 3DEP database."""
import tempfile
import uuid
from pathlib import Path
from typing import Dict, List, Tuple, Union

import async_retriever as ar
import cytoolz as tlz
import numpy as np
import pygeoutils as geoutils
import pyproj
import rasterio as rio
import rasterio.warp as rio_warp
import xarray as xr
from pydantic import BaseModel, validator
from pygeoogc import WMS, InvalidInputType, InvalidInputValue, ServiceURL, utils
from shapely.geometry import MultiPolygon, Polygon

try:
    import richdem as rd
except ImportError:
    rd = None

from .exceptions import MissingDependency

DEF_CRS = "epsg:4326"
LAYERS = [
    "DEM",
    "Hillshade Gray",
    "Aspect Degrees",
    "Aspect Map",
    "GreyHillshade_elevationFill",
    "Hillshade Multidirectional",
    "Slope Map",
    "Slope Degrees",
    "Hillshade Elevation Tinted",
    "Height Ellipsoidal",
    "Contour 25",
    "Contour Smoothed 25",
]
__all__ = ["get_map", "elevation_bygrid", "elevation_bycoords", "deg2mpm"]


def get_map(
    layers: Union[str, List[str]],
    geometry: Union[Polygon, MultiPolygon, Tuple[float, float, float, float]],
    resolution: float,
    geo_crs: str = DEF_CRS,
    crs: str = DEF_CRS,
) -> Union[xr.DataArray, xr.Dataset]:
    """Access to `3DEP <https://www.usgs.gov/core-science-systems/ngp/3dep>`__ service.

    The 3DEP service has multi-resolution sources so depending on the user
    provided resolution the data is resampled on server-side based
    on all the available data sources. The following layers are available:

    - ``DEM``
    - ``Hillshade Gray``
    - ``Aspect Degrees``
    - ``Aspect Map``
    - ``GreyHillshade_elevationFill``
    - ``Hillshade Multidirectional``
    - ``Slope Map``
    - ``Slope Degrees``
    - ``Hillshade Elevation Tinted``
    - ``Height Ellipsoidal``
    - ``Contour 25``
    - ``Contour Smoothed 25``

    Parameters
    ----------
    layers : str or list
        A valid 3DEP layer or a list of them.
    geometry : Polygon, MultiPolygon, or tuple
        A shapely Polygon or a bounding box ``(west, south, east, north)``.
    resolution : float
        The data resolution in meters. The width and height of the output are computed in
        pixels
        based on the geometry bounds and the given resolution.
    geo_crs : str, optional
        The spatial reference system of the input geometry, defaults to
        ``EPSG:4326``.
    crs : str, optional
        The spatial reference system to be used for requesting the data, defaults to
        ``EPSG:4326``. Valis values are ``epsg:4326``, ``epsg:3576``, ``epsg:3571``,
        ``epsg:3575``, ``epsg:3857``, ``epsg:3572``, ``crs:84``, ``epsg:3573``,
        and ``epsg:3574``.

    Returns
    -------
    dict
        A dict where the keys are the layer name and values are the returned response
        from the WMS service as bytes. You can use ``utils.create_dataset`` function
        to convert the responses to ``xarray.Dataset``.
    """
    _layers = layers.copy() if isinstance(layers, list) else [layers]
    invalid = [lyr for lyr in _layers if lyr not in LAYERS]
    if len(invalid) > 0:
        raise InvalidInputValue("layers", LAYERS)

    if "DEM" in _layers:
        _layers[_layers.index("DEM")] = "None"

    _layers = [f"3DEPElevation:{lyr}" for lyr in _layers]

    _geometry = geoutils.geo2polygon(geometry, geo_crs, crs)
    wms = WMS(ServiceURL().wms.nm_3dep, layers=_layers, outformat="image/tiff", crs=crs)
    r_dict = wms.getmap_bybox(_geometry.bounds, resolution, box_crs=crs)

    ds = geoutils.gtiff2xarray(r_dict, _geometry, crs)

    valid_layers = wms.get_validlayers()
    rename = {lyr: lyr.split(":")[-1].replace(" ", "_").lower() for lyr in valid_layers}
    rename.update({"3DEPElevation:None": "elevation"})

    if isinstance(ds, xr.DataArray):
        ds.name = rename[ds.name]
    else:
        ds = ds.rename({n: rename[n] for n in ds.keys()})

    return ds


def elevation_bygrid(
    xcoords: List[float],
    ycoords: List[float],
    crs: str,
    resolution: float,
    depression_filling: bool = False,
) -> xr.DataArray:
    """Get elevation from DEM data for a grid.

    This function is intended for getting elevations for a gridded dataset.

    Parameters
    ----------
    xcoords : list
        List x-coordinates of a a grid.
    ycoords : list
        List of y-coordinates of a grid.
    crs : str
        The spatial reference system of the input grid, defaults to ``EPSG:4326``.
    resolution : float
        The accuracy of the output, defaults to 10 m which is the highest
        available resolution that covers CONUS. Note that higher resolution
        increases computation time so chose this value with caution.
    depression_filling : bool, optional
        Fill depressions before sampling using
        `RichDEM <https://richdem.readthedocs.io/en/latest/>`__ package, defaults to False.

    Returns
    -------
    xarray.DataArray
        An data array with name elevation and the given dim names.
    """
    bbox = (min(xcoords), min(ycoords), max(xcoords), max(ycoords))

    ratio_min = 0.01
    ratio_x = abs((bbox[2] - bbox[0]) / bbox[0])
    ratio_y = abs((bbox[3] - bbox[1]) / bbox[1])
    if (ratio_x < ratio_min) or (ratio_y < ratio_min):
        rad = ratio_min * abs(bbox[0])
        bbox = (bbox[0] - rad, bbox[1] - rad, bbox[2] + rad, bbox[3] + rad)

    wms = WMS(
        ServiceURL().wms.nm_3dep, layers="3DEPElevation:None", outformat="image/tiff", crs=DEF_CRS
    )
    r_dict = wms.getmap_bybox(bbox, resolution, box_crs=crs)
    dem = _reproject_gtiff(r_dict, crs)

    if depression_filling:
        dem = fill_depressions(dem)

    return dem.interp(x=xcoords, y=ycoords)


def _reproject_gtiff(
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
        Conditioned DEM.
    """
    if rd is None:
        raise MissingDependency

    rda = rd.rdarray(dem, no_data=dem.nodatavals[0])
    rda.projection = dem.crs
    rda.geotransform = dem.transform
    rda = rd.FillDepressions(rda, epsilon=False)
    dem.data = rd.ResolveFlats(rda)
    return dem


class ElevationByCoords(BaseModel):
    """Elevation model class.

    Parameters
    ----------
    coords : list of tuple
        List of coordinates.
    crs : str, optional
        Coordinate reference system of the input coordinates, defaults to ``EPSG:4326``.
    source : str, optional
        Elevation source, defaults to ``tnm``. Valid sources are: ``tnm``, ``airmap``.
    """

    crs: str = DEF_CRS
    coords: List[Tuple[float, float]]
    source: str = "tnm"

    @validator("crs")
    def _valid_crs(cls, v):
        try:
            return pyproj.CRS(v)
        except pyproj.exceptions.CRSError as ex:
            raise InvalidInputType("crs", "a valid CRS") from ex

    @validator("coords")
    def _validate_coords(cls, v, values):
        return utils.match_crs(v, values["crs"], DEF_CRS)

    @validator("source")
    def _validate_source(cls, v):
        valid_sources = ["tnm", "airmap"]
        if v not in valid_sources:
            raise InvalidInputValue("source", valid_sources)
        return v

    def airmap(self) -> List[Tuple[float, float]]:
        """Return list of elevations in meters."""
        coords_chunks = tlz.partition_all(100, self.coords)
        headers = {"Content-Type": "application/json", "charset": "utf-8"}
        urls, kwds = zip(
            *(
                (
                    ServiceURL().restful.airmap,
                    {
                        "params": {"points": ",".join(f"{lat},{lon}" for lon, lat in chunk)},
                        "headers": headers,
                    },
                )
                for chunk in coords_chunks
            )
        )
        elevations = list(tlz.pluck("data", ar.retrieve(urls, "json", kwds)))
        return list(tlz.concat(elevations))

    def tnm(self) -> List[Tuple[float, float]]:
        """Return list of elevations in meters."""
        urls, kwds = zip(
            *(
                (
                    ServiceURL().restful.nm_pqs,
                    {
                        "params": {
                            "units": "meters",
                            "output": "json",
                            "y": f"{lat:.5f}",
                            "x": f"{lon:.5f}",
                        },
                    },
                )
                for lon, lat in self.coords
            )
        )
        resp = ar.retrieve(urls, "json", kwds, max_workers=200)
        return [
            r["USGS_Elevation_Point_Query_Service"]["Elevation_Query"]["Elevation"] for r in resp
        ]


def elevation_bycoords(
    coords: List[Tuple[float, float]], crs: Union[str, pyproj.CRS] = DEF_CRS, source: str = "tnm"
) -> List[float]:
    """Get elevation from Airmap at 1-arc resolution (~30 m) for a list of coordinates.

    Parameters
    ----------
    coords : list of tuples
        Coordinates of target location as list of tuples ``[(x, y), ...]``.
    crs : str or pyproj.CRS, optional
        Spatial reference (CRS) of coords, defaults to ``EPSG:4326``.
    source : str, optional
        Data source to be used, default to ``tnm``. Supported sources are
        ``airmap`` (30 m resolution) and ``tnm`` (using The National Map's Bulk Point
        Query Service). The ``tnm`` source is more accurate since it uses the highest
        available resolution DEM automatically but it is limited to the US.

    Returns
    -------
    list of float
        Elevation in meter.
    """
    _crs = crs.to_string() if isinstance(crs, pyproj.CRS) else crs
    service = ElevationByCoords(crs=_crs, coords=coords, source=source)
    return getattr(service, source)()


def deg2mpm(da: xr.DataArray) -> xr.DataArray:
    """Convert slope from degree to meter/meter."""
    attrs = da.attrs
    da = np.tan(np.deg2rad(da))
    da.attrs = attrs
    da.name = "slope"
    da.attrs["units"] = "meters/meters"
    return da
