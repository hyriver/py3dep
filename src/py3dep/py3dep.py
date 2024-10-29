"""Get data from 3DEP database."""

# pyright: reportGeneralTypeIssues=false
from __future__ import annotations

import contextlib
import itertools
from collections.abc import Sequence
from typing import TYPE_CHECKING, Any, Literal, Union, cast, overload

import geopandas as gpd
import numpy as np
import pandas as pd
import pyproj
import rasterio
import rioxarray as rxr
import xarray as xr
from rasterio import RasterioIOError
from shapely import LineString, MultiLineString, MultiPolygon, Polygon, ops
from shapely import box as shapely_box

import async_retriever as ar
import pygeoutils as geoutils
from py3dep import utils
from py3dep.exceptions import (
    InputTypeError,
    InputValueError,
    MissingCRSError,
    ServiceUnavailableError,
)
from pygeoogc import WMS, ArcGISRESTful, ServiceURL
from pygeoogc import utils as ogc_utils
from pygeoogc.exceptions import ZeroMatchedError

if TYPE_CHECKING:
    from pathlib import Path

    CRSTYPE = Union[int, str, pyproj.CRS]

MAX_PIXELS = 8000000
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
__all__ = [
    "get_map",
    "elevation_bygrid",
    "elevation_bycoords",
    "elevation_profile",
    "check_3dep_availability",
    "query_3dep_sources",
    "static_3dep_dem",
    "get_dem",
    "get_dem_vrt",
    "add_elevation",
]


@overload
def get_map(
    layers: str,
    geometry: Polygon | MultiPolygon | tuple[float, float, float, float],
    resolution: int,
    geo_crs: CRSTYPE = ...,
    crs: CRSTYPE = ...,
) -> xr.DataArray: ...


@overload
def get_map(
    layers: list[str],
    geometry: Polygon | MultiPolygon | tuple[float, float, float, float],
    resolution: int,
    geo_crs: CRSTYPE = ...,
    crs: CRSTYPE = ...,
) -> xr.Dataset: ...


def get_map(
    layers: str | list[str],
    geometry: Polygon | MultiPolygon | tuple[float, float, float, float],
    resolution: int,
    geo_crs: CRSTYPE = 4326,
    crs: CRSTYPE = 4326,
) -> xr.Dataset | xr.DataArray:
    """Access dynamic layer of `3DEP <https://www.usgs.gov/core-science-systems/ngp/3dep>`__.

    The 3DEP service has multi-resolution sources, so depending on the user
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
    layers : str or list of str
        A valid 3DEP layer or a list of them.
    geometry : Polygon, MultiPolygon, or tuple
        A shapely Polygon or a bounding box of the form ``(west, south, east, north)``.
    resolution : int
        The target resolution in meters. The width and height of the output are computed in
        pixels based on the geometry bounds and the given resolution.
    geo_crs : str, int, or pyproj.CRS, optional
        The spatial reference system of the input geometry, defaults to ``EPSG:4326``.
    crs : str, int, or pyproj.CRS, optional
        The spatial reference system to be used for requesting the data, defaults to
        ``EPSG:4326``. Valid values are ``EPSG:4326``, ``EPSG:3576``, ``EPSG:3571``,
        ``EPSG:3575``, ``EPSG:3857``, ``EPSG:3572``, ``CRS:84``, ``EPSG:3573``,
        and ``EPSG:3574``.

    Returns
    -------
    xarray.DataArray or xarray.Dataset
        The requested topographic data as an ``xarray.DataArray`` or ``xarray.Dataset``.
    """
    _layers = list(layers) if isinstance(layers, (list, tuple)) else [layers]
    invalid = [lyr for lyr in _layers if lyr not in LAYERS]
    if invalid:
        raise InputValueError("layers", LAYERS, ", ".join(invalid))

    if "DEM" in _layers:
        _layers[_layers.index("DEM")] = "None"

    _layers = [f"3DEPElevation:{lyr}" for lyr in _layers]

    wms_url = ServiceURL().wms.nm_3dep

    valid_crs = ogc_utils.valid_wms_crs(wms_url)

    if len(valid_crs) == 0:
        raise ServiceUnavailableError(wms_url)

    if ogc_utils.validate_crs(crs).lower() not in valid_crs:
        raise InputValueError("crs", valid_crs)

    _geometry = geoutils.geo2polygon(geometry, geo_crs, crs)
    wms = WMS(wms_url, layers=_layers, outformat="image/tiff", crs=crs, validation=False)
    r_dict = wms.getmap_bybox(_geometry.bounds, resolution, box_crs=crs, max_px=MAX_PIXELS)

    try:
        ds = geoutils.gtiff2xarray(r_dict, _geometry, crs)
    except RasterioIOError as ex:
        raise ServiceUnavailableError(wms_url) from ex
    valid_layers = wms.get_validlayers()
    return utils.rename_layers(ds, list(valid_layers))


def static_3dep_dem(
    geometry: Polygon | MultiPolygon | tuple[float, float, float, float],
    crs: CRSTYPE,
    resolution: int = 10,
) -> xr.DataArray:
    """Get DEM data at specific resolution from 3DEP.

    Notes
    -----
    In contrast to ``get_map`` function, this function only gets DEM data at
    specific resolution, namely 10 m, 30 m, and 60 m. However, this function
    is faster. This function is intended for cases where only need DEM at a
    specific resolution is required and for the other requests ``get_map``
    should be used.

    Parameters
    ----------
    geometry : Polygon, MultiPolygon, or tuple of length 4
        Geometry to get DEM within. It can be a polygon or a boundong box
        of form (xmin, ymin, xmax, ymax).
    crs : int, str, of pyproj.CRS
        CRS of the input geometry.
    resolution : int, optional
        Target DEM source resolution in meters, defaults to 10 m which is the highest
        resolution available over the US. Available options are 10, 30, and 60.

    Returns
    -------
    xarray.DataArray
        The request DEM at the specified resolution.
    """
    base_url = "https://prd-tnm.s3.amazonaws.com/StagedProducts/Elevation"
    url = {
        10: f"{base_url}/13/TIFF/USGS_Seamless_DEM_13.vrt",
        30: f"{base_url}/1/TIFF/USGS_Seamless_DEM_1.vrt",
        60: f"{base_url}/2/TIFF/USGS_Seamless_DEM_2.vrt",
    }
    if resolution not in url:
        raise InputValueError("resolution", list(url))

    dem = rxr.open_rasterio(url[resolution])
    dem = cast("xr.DataArray", dem)
    if "band" in dem.dims:
        dem = dem.squeeze("band")
    poly = geoutils.geo2polygon(geometry, crs, dem.rio.crs)
    dem = dem.rio.clip_box(*poly.bounds)
    if isinstance(geometry, (Polygon, MultiPolygon)):
        dem = dem.rio.clip([poly])
    dem = dem.where(dem > dem.rio.nodata, drop=False)
    dem = dem.rio.write_nodata(np.nan)
    dem.attrs.update({"units": "meters", "vertical_datum": "NAVD88", "vertical_resolution": 0.001})
    dem.name = "elevation"
    return dem


def get_dem(
    geometry: Polygon | MultiPolygon | tuple[float, float, float, float],
    resolution: int,
    crs: CRSTYPE = 4326,
) -> xr.DataArray:
    """Get DEM data at any resolution from 3DEP.

    Notes
    -----
    This function is a wrapper of ``static_3dep_dem`` and ``get_map`` functions.
    Since ``static_3dep_dem`` is much faster, if the requested resolution is 10 m,
    30 m, or 60 m, ``static_3dep_dem`` will be used. Otherwise, ``get_map``
    will be used.

    Parameters
    ----------
    geometry : Polygon, MultiPolygon, or tuple of length 4
        Geometry to get DEM within. It can be a polygon or a boundong box
        of form (xmin, ymin, xmax, ymax).
    resolution : int
        Target DEM source resolution in meters.
    crs : str, int, or pyproj.CRS, optional
        The spatial reference system of the input geometry, defaults to ``EPSG:4326``.

    Returns
    -------
    xarray.DataArray
        DEM at the specified resolution in meters and 4326 CRS.
    """
    if np.isclose(resolution, (10, 30, 60)).any():
        dem = static_3dep_dem(geometry, crs, resolution)
    else:
        dem = get_map("DEM", geometry, resolution, crs)
    dem = dem.astype("f4")
    dem.attrs.update({"units": "m", "vertical_datum": "NAVD88"})
    dem.name = "elevation"
    return dem


def add_elevation(
    ds: xr.DataArray | xr.Dataset,
    resolution: int | None = None,
    x_dim: str = "x",
    y_dim: str = "y",
    mask: xr.DataArray | None = None,
) -> xr.Dataset:
    """Add elevation data to a dataset as a new variable.

    Parameters
    ----------
    ds : xarray.DataArray or xarray.Dataset
        The dataset to add elevation data to. It must contain
        CRS information.
    resolution : float, optional
        Target DEM source resolution in meters, defaults ``None``, i.e.,
        the resolution of the input ``ds`` will be used.
    x_dim : str, optional
        Name of the x-coordinate dimension in ``ds``, defaults to ``x``.
    y_dim : str, optional
        Name of the y-coordinate dimension in ``ds``, defaults to ``y``.
    mask : xarray.DataArray, optional
        A mask to apply to the elevation data, defaults to ``None``.

    Returns
    -------
    xarray.Dataset
        The dataset with ``elevation`` variable added.
    """
    if not isinstance(ds, (xr.DataArray, xr.Dataset)):
        raise InputTypeError("ds", "xarray.DataArray or xarray.Dataset")

    ds_crs = ds.rio.crs
    if ds_crs is None:
        raise MissingCRSError

    if isinstance(ds, xr.DataArray):
        name = ds.name or "data"
        ds = ds.to_dataset(name=name, promote_attrs=True)
    else:
        ds = ds.copy()
    ds = ds.rio.set_spatial_dims(x_dim=x_dim, y_dim=y_dim)
    if ds_crs.is_projected:
        ds_proj = ds
        crs_proj = ds_crs
    else:
        ds_proj = ds.rio.reproject(5070)
        crs_proj = 5070
    ds_bounds = ds_proj.rio.bounds()
    if resolution is None:
        resolution = int(ds_proj.rio.resolution()[0])
    bounds = shapely_box(*ds_bounds).buffer(3 * resolution, join_style=2, cap_style=2)
    bounds = geoutils.geometry_reproject(bounds, crs_proj, 4326)

    dims_order = (y_dim, x_dim)
    elev = (
        get_dem(bounds.bounds, resolution).rename({"x": x_dim, "y": y_dim}).transpose(*dims_order)
    )
    if ds_crs != elev.rio.crs:
        elev = elev.rio.reproject(ds_crs)
    elev = geoutils.xarray_geomask(elev, ds.rio.bounds(), ds_crs)
    ds["elevation"] = (dims_order, elev.rio.reproject_match(ds, resampling=1).values)
    ds["elevation"] = ds["elevation"].rio.write_crs(
        ds.rio.crs, grid_mapping_name=ds.rio.grid_mapping
    )
    if mask is not None:
        ds["elevation"] = ds["elevation"].where(mask)
    ds = cast("xr.Dataset", ds)
    return ds


def get_dem_vrt(
    bbox: tuple[float, float, float, float],
    resolution: int,
    vrt_path: Path | str,
    tiff_dir: Path | str = "cache",
    crs: CRSTYPE = 4326,
) -> None:
    """Get DEM data at any resolution from 3DEP and save it as a VRT file.

    Parameters
    ----------
    bbox : tuple of length 4
        The boundong box of form (xmin, ymin, xmax, ymax).
    resolution : int
        Target DEM source resolution in meters.
    vrt_path : str or pathlib.Path
        Path to the output VRT file.
    tiff_dir : str or pathlib.Path, optional
        Path to the directory to save the downloaded TIFF file, defaults
        to ``./cache``.
    crs : str, int, or pyproj.CRS, optional
        The spatial reference system of ``bbox``, defaults to ``EPSG:4326``.
    """
    wms_url = ServiceURL().wms.nm_3dep
    bounds = geoutils.geometry_reproject(bbox, crs, 4326)
    fname = WMS(
        wms_url, layers="3DEPElevation:None", outformat="image/tiff", crs=4326, validation=False
    ).getmap_bybox(bounds, resolution, max_px=MAX_PIXELS, tiff_dir=tiff_dir)
    return geoutils.gtiff2vrt(fname, vrt_path)


def elevation_bygrid(
    xcoords: Sequence[float],
    ycoords: Sequence[float],
    crs: CRSTYPE,
    resolution: int,
    depression_filling: bool = False,
) -> xr.DataArray:
    """Get elevation from DEM data for a grid.

    This function is intended for getting elevations for a gridded dataset.

    Parameters
    ----------
    xcoords : list
        List of x-coordinates of a grid.
    ycoords : list
        List of y-coordinates of a grid.
    crs : str, int, or pyproj.CRS or pyproj.CRS
        The spatial reference system of the input grid,
        defaults to ``EPSG:4326``.
    resolution : int
        The accuracy of the output, defaults to 10 m which is the highest
        available resolution that covers CONUS. Note that higher resolution
        increases computation time so chose this value with caution.
    depression_filling : bool, optional
        Fill depressions before sampling using
        `Wang and Liu (2006) <https://doi.org/10.1080/13658810500433453>`__
        method, defaults to ``False``.

    Returns
    -------
    xarray.DataArray
        Elevations of the input coordinates as a ``xarray.DataArray``.
    """
    pts_crs = ogc_utils.validate_crs(crs)
    bbox = (
        gpd.GeoSeries(gpd.points_from_xy(*zip(*itertools.product(xcoords, ycoords)), crs=pts_crs))
        .to_crs(5070)
        .buffer(2 * resolution)
        .total_bounds
    )

    dem = get_dem(tuple(bbox), resolution, 5070)
    dem = dem.rio.reproject(pts_crs)

    if depression_filling:
        dem = utils.fill_depressions(dem)

    return dem.interp(x=list(xcoords), y=list(ycoords))


class ElevationByCoords:
    """Elevation model class.

    Parameters
    ----------
    crs : str, int, or pyproj.CRS, optional
        Coordinate reference system of the input coordinates, defaults to ``EPSG:4326``.
    coords : list of tuple
        List of coordinates.
    source : str, optional
        Elevation source, defaults to ``tep``. Valid sources are: ``tnm``, and ``tep``.
    """

    def __init__(
        self,
        coords: list[tuple[float, float]] | tuple[float, float],
        crs: CRSTYPE = 4326,
        source: str = "tep",
    ) -> None:
        self.coords = geoutils.coords_list(coords)
        self.crs = ogc_utils.validate_crs(crs)
        self.coords_gs = gpd.GeoSeries(gpd.points_from_xy(*zip(*self.coords)), crs=self.crs)
        valid_sources = ("tnm", "tep")
        self.source = source
        if self.source not in valid_sources:
            raise InputValueError("source", valid_sources)

    @property
    def values(self) -> list[float]:
        """Return list of elevations in meters."""
        if self.source == "tep":
            return self.tep()
        return self.tnm()

    def tnm(self) -> list[float]:
        """Return list of elevations in meters."""
        pts = self.coords_gs.to_crs(4326)
        kwds = [
            {"params": {"units": "Meters", "x": f"{x:.6f}", "y": f"{y:.6f}"}}
            for x, y in zip(pts.x, pts.y)
        ]
        resp = ar.retrieve_json([ServiceURL().restful.nm_pqs] * len(kwds), kwds, max_workers=5)
        resp = cast("list[dict[str, Any]]", resp)
        return [float(r["value"]) for r in resp]

    def tep(self) -> list[float]:
        """Get elevation from 10-m resolution 3DEP."""
        url = "/".join(
            (
                "https://prd-tnm.s3.amazonaws.com/StagedProducts/Elevation",
                "13/TIFF/USGS_Seamless_DEM_13.vrt",
            )
        )
        with rasterio.open(url) as src:
            points_proj = self.coords_gs.to_crs(src.crs)
            elev = np.array(list(src.sample(zip(points_proj.x, points_proj.y)))).ravel()
        return elev.tolist()


@overload
def elevation_bycoords(
    coords: tuple[float, float], crs: CRSTYPE = ..., source: Literal["tep", "tnm"] = ...
) -> float: ...


@overload
def elevation_bycoords(
    coords: list[tuple[float, float]],
    crs: CRSTYPE = ...,
    source: Literal["tep", "tnm"] = ...,
) -> list[float]: ...


def elevation_bycoords(
    coords: tuple[float, float] | list[tuple[float, float]],
    crs: CRSTYPE = 4326,
    source: Literal["tep", "tnm"] = "tep",
) -> float | list[float]:
    """Get elevation for a list of coordinates.

    Parameters
    ----------
    coords : tuple or list of tuple
        Coordinates of target location(s), e.g., ``[(x, y), ...]``.
    crs : str, int, or pyproj.CRS or pyproj.CRS, optional
        Spatial reference (CRS) of coords, defaults to ``EPSG:4326``.
    source : str, optional
        Data source to be used, default to ``tep``. Supported sources are
        ``tnm`` (using The National Map's Bulk Point
        Query Service with 10 m resolution) and ``tep`` (using 3DEP's static DEM VRTs
        at 10 m resolution). The ``tnm`` and ``tep`` sources are more accurate since they
        use the 1/3 arc-second DEM layer from 3DEP service but it is limited to the US.
        Note that ``tnm`` is bit unstable. It's recommended to use ``tep`` unless 10-m
        resolution accuracy is not necessary.

    Returns
    -------
    float or list of float
        Elevation in meter.
    """
    _crs = crs.to_string() if isinstance(crs, pyproj.CRS) else crs
    service = ElevationByCoords(crs=_crs, coords=coords, source=source)
    if len(service.coords) == 1:
        return service.values[0]
    return service.values


def elevation_profile(
    lines: LineString | MultiLineString,
    spacing: float,
    crs: CRSTYPE = 4326,
) -> xr.DataArray:
    """Get the elevation profile along a line at a given uniform spacing.

    .. note::

        This function converts the line to a spline and then calculates the elevation
        along the spline at a given uniform spacing using 10-m resolution DEM from 3DEP.

    Parameters
    ----------
    lines : LineString or MultiLineString
        Line segment(s) to be profiled. If its type is ``MultiLineString``,
        it will be converted to a single ``LineString`` and if this operation
        fails, an ``InputTypeError`` will be raised.
    spacing : float
        Spacing between the sample points along the line in meters.
    crs : str, int, or pyproj.CRS, optional
        Spatial reference System (CRS) of ``lines``, defaults to ``EPSG:4326``.

    Returns
    -------
    xarray.DataArray
        Elevation profile with dimension ``z`` and three coordinates: ``x``, ``y``,
        and ``distance``. The ``distance`` coordinate is the distance from the start
        of the line in meters.
    """
    if not isinstance(lines, (LineString, MultiLineString)):
        raise InputTypeError("lines", "LineString or MultiLineString")

    if isinstance(lines, MultiLineString):
        path = ops.linemerge(lines)
        if not isinstance(path, LineString):
            raise InputTypeError("lines", "MultiLineString that can be merged to LineString")
    else:
        path = lines

    crs_prj = 5070
    geom_proj = geoutils.geometry_reproject(path, crs, crs_prj)
    npts = int(np.ceil(geom_proj.length / spacing))
    geom_proj = geoutils.smooth_linestring(geom_proj, 0.1, npts)
    elev_list = elevation_bycoords(list(zip(*geom_proj.xy)), crs=crs_prj, source="tep")
    elevation = xr.DataArray(
        elev_list,
        dims="z",
        coords={"z": range(len(elev_list))},
        attrs={"source": "10-m DEM from 3DEP"},
    )
    geom = geoutils.geometry_reproject(geom_proj, crs_prj, crs)
    x, y = geom.xy
    elevation["x"], elevation["y"] = ("z", list(x)), ("z", list(y))
    distance = np.hypot(np.diff(x), np.diff(y)).cumsum()
    distance = np.insert(distance, 0, 0)
    elevation["distance"] = ("z", distance)
    elevation["distance"].attrs = {"units": "m", "long_name": "Distance from start"}
    return elevation


def check_3dep_availability(
    bbox: tuple[float, float, float, float], crs: CRSTYPE = 4326
) -> dict[str, bool | str]:
    """Query 3DEP's resolution availability within a bounding box.

    This function checks availability of 3DEP's at the following resolutions:
    1 m, 3 m, 5 m, 10 m, 30 m, 60 m, and topobathy (integrated topobathymetry).

    Parameters
    ----------
    bbox : tuple
        Bounding box as tuple of ``(min_x, min_y, max_x, max_y)``.
    crs : str, int, or pyproj.CRS or pyproj.CRS, optional
        Spatial reference (CRS) of ``bbox``, defaults to ``EPSG:4326``.

    Returns
    -------
    dict
        ``True`` if bbox intersects 3DEP elevation for each available resolution.
        Keys are the supported resolutions and values are their availability.
        If the query fails due to any reason, the value will be ``Failed``.
        If necessary, you can try again later until there is no ``Failed`` value.

    Examples
    --------
    >>> import py3dep
    >>> bbox = (-69.77, 45.07, -69.31, 45.45)
    >>> py3dep.check_3dep_availability(bbox)
    {'1m': True, '3m': False, '5m': False, '10m': True, '30m': True, '60m': False, 'topobathy': False}
    """
    if not isinstance(bbox, Sequence) or len(bbox) != 4:
        raise InputTypeError("bbox", "a tuple of length 4")

    res_layers = {
        "1m": 18,
        "3m": 19,
        "5m": 20,
        "10m": 21,
        "30m": 22,
        "60m": 23,
        "topobathy": 30,
    }
    base_url = ServiceURL().restful.nm_3dep_index
    urls = [f"{base_url}/{lyr}/query" for lyr in res_layers.values()]

    geom_query = ogc_utils.esri_query(bbox, crs, 4326)
    payload = {
        **geom_query,
        "spatialRel": "esriSpatialRelIntersects",
        "returnGeometry": "false",
        "returnIdsOnly": "false",
        "returnCountOnly": "true",
        "f": "json",
    }

    resps = ar.retrieve_json(
        urls,
        [{"params": payload}] * len(urls),
        max_workers=len(urls),
        raise_status=False,
    )
    resps = cast("list[dict[str, Any]]", resps)

    def _check(r: dict[str, Any] | None) -> bool | str:
        if r is None or "error" in r:
            return "Failed"
        return bool(r.get("count"))

    avail = {res: _check(r) for res, r in zip(res_layers, resps)}
    failed = [res for res, r in avail.items() if r == "Failed"]
    if failed:
        [
            ar.delete_url_cache(f"{base_url}/{res_layers[res]}/query", params=payload)
            for res in failed
        ]
    return avail


def query_3dep_sources(
    bbox: tuple[float, float, float, float],
    crs: CRSTYPE = 4326,
    res: str | list[str] | None = None,
) -> gpd.GeoDataFrame:
    """Query 3DEP's data sources within a bounding box.

    This function queries the availability of the underlying data that 3DEP uses
    at the following resolutions:
    1 m, 3 m, 5 m, 10 m, 30 m, 60 m, and topobathy (integrated topobathymetry).

    Parameters
    ----------
    bbox : tuple
        Bounding box as tuple of ``(min_x, min_y, max_x, max_y)``.
    crs : str, int, or pyproj.CRS or pyproj.CRS, optional
        Spatial reference (CRS) of bbox, defaults to ``EPSG:4326``.
    res : str, list of str, optional
        Resolution to query, defaults to ``None``, i.e., all resolutions.
        Available resolutions are: ``1m``, ``3m``, ``5m``, ``10m``, ``30m``,
        ``60m``, and ``topobathy``.

    Returns
    -------
    geopandas.GeoDataFrame
        Polygon(s) representing the 3DEP data sources at each resolution.
        Resolutions are given in the ``dem_res`` column.

    Examples
    --------
    >>> import py3dep
    >>> bbox = (-69.77, 45.07, -69.31, 45.45)
    >>> src = py3dep.query_3dep_sources(bbox)
    >>> src.groupby("dem_res")["OBJECTID"].count().to_dict()
    {'10m': 16, '1m': 4, '30m': 8}
    >>> src = py3dep.query_3dep_sources(bbox, res="1m")
    >>> src.groupby("dem_res")["OBJECTID"].count().to_dict()
    {'1m': 4}
    """
    if not isinstance(bbox, Sequence) or len(bbox) != 4:
        raise InputTypeError("bbox", "a tuple of length 4")

    res_layers = {
        "1m": 18,
        "3m": 19,
        "5m": 20,
        "10m": 21,
        "30m": 22,
        "60m": 23,
        "topobathy": 30,
    }
    if res is None:
        layers = res_layers
    elif isinstance(res, str) and res in res_layers:
        layers = {res: res_layers[res]}
    elif isinstance(res, (list, tuple)) and all(r in res_layers for r in res):
        layers = {r: res_layers[r] for r in res}
    else:
        raise InputValueError("res", list(res_layers))

    def _check(lyr: int) -> gpd.GeoDataFrame | None:
        client = ArcGISRESTful(ServiceURL().restful.nm_3dep_index, lyr, outformat="json")
        with contextlib.suppress(ZeroMatchedError):
            oids = client.oids_bygeom(bbox, crs)
            return geoutils.json2geodf(client.get_features(oids))
        return None

    return (  # pyright: ignore[reportReturnType]
        gpd.GeoDataFrame(  # pyright: ignore[reportCallIssue]
            pd.concat({res: _check(lyr) for res, lyr in layers.items()}).reset_index(  # pyright: ignore[reportCallIssue,reportArgumentType]
                level=1, drop=True
            ),
            crs=4326,
        )
        .reset_index()
        .rename(columns={"index": "dem_res"})
    )
