"""Get data from 3DEP database."""
import contextlib
from typing import Dict, List, Sequence, Tuple, Union

import async_retriever as ar
import cytoolz as tlz
import geopandas as gpd
import pygeoutils as geoutils
import pyproj
import xarray as xr
from pydantic import BaseModel, validator
from pygeoogc import (
    WMS,
    ArcGISRESTful,
    InvalidInputType,
    InvalidInputValue,
    ServiceError,
    ServiceURL,
)
from pygeoogc import utils as ogc_utils
from shapely.geometry import MultiPolygon, Polygon

from . import utils

DEF_CRS = "epsg:4326"
EXPIRE = -1
LAYERS = [
    "DEM",
    "Hillshade Gray",
    "Aspect Degrees",
    "Aspect Map",
    "GreyHillshade_elevationFill",
    "Hillshade Multidirectional",
    "Slope Map",
    "Hillshade Elevation Tinted",
    "Height Ellipsoidal",
    "Contour 25",
    "Contour Smoothed 25",
]
__all__ = ["get_map", "elevation_bygrid", "elevation_bycoords", "check_3dep_availability"]


def get_map(
    layers: Union[str, Sequence[str]],
    geometry: Union[Polygon, MultiPolygon, Tuple[float, float, float, float]],
    resolution: float,
    geo_crs: str = DEF_CRS,
    crs: str = DEF_CRS,
    expire_after: int = EXPIRE,
    disable_caching: bool = False,
) -> Union[xr.Dataset, xr.DataArray]:
    """Access to `3DEP <https://www.usgs.gov/core-science-systems/ngp/3dep>`__ service.

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
    resolution : float
        The target resolution in meters. The width and height of the output are computed in
        pixels based on the geometry bounds and the given resolution.
    geo_crs : str, optional
        The spatial reference system of the input geometry, defaults to ``EPSG:4326``.
    crs : str, optional
        The spatial reference system to be used for requesting the data, defaults to
        ``EPSG:4326``. Valid values are ``EPSG:4326``, ``EPSG:3576``, ``EPSG:3571``,
        ``EPSG:3575``, ``EPSG:3857``, ``EPSG:3572``, ``CRS:84``, ``EPSG:3573``,
        and ``EPSG:3574``.
    expire_after : int, optional
        Expiration time for response caching in seconds, defaults to -1 (never expire).
    disable_caching : bool, optional
        If ``True``, disable caching requests, defaults to False.

    Returns
    -------
    xarray.DataArray or xarray.Dataset
        The requested topographic data as an ``xarray.DataArray`` or ``xarray.Dataset``.
    """
    _layers = list(layers) if isinstance(layers, (tuple, list)) else [layers]
    invalid = [lyr for lyr in _layers if lyr not in LAYERS]
    if len(invalid) > 0:
        raise InvalidInputValue("layers", LAYERS)

    if "DEM" in _layers:
        _layers[_layers.index("DEM")] = "None"

    _layers = [f"3DEPElevation:{lyr}" for lyr in _layers]

    _geometry = geoutils.geo2polygon(geometry, geo_crs, crs)
    wms = WMS(
        ServiceURL().wms.nm_3dep,
        layers=_layers,
        outformat="image/tiff",
        crs=crs,
        expire_after=expire_after,
        disable_caching=disable_caching,
    )
    r_dict = wms.getmap_bybox(_geometry.bounds, resolution, box_crs=crs)
    if isinstance(geometry, (Polygon, MultiPolygon)):
        ds: xr.Dataset = geoutils.gtiff2xarray(r_dict, _geometry, crs)
    else:
        ds = geoutils.gtiff2xarray(r_dict)
    valid_layers = wms.get_validlayers()
    return utils.rename_layers(ds, list(valid_layers))


def elevation_bygrid(
    xcoords: Sequence[float],
    ycoords: Sequence[float],
    crs: str,
    resolution: float,
    depression_filling: bool = False,
    expire_after: int = EXPIRE,
    disable_caching: bool = False,
) -> xr.DataArray:
    """Get elevation from DEM data for a grid.

    This function is intended for getting elevations for a gridded dataset.

    Parameters
    ----------
    xcoords : list
        List of x-coordinates of a grid.
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
    expire_after : int, optional
        Expiration time for response caching in seconds, defaults to -1 (never expire).
    disable_caching : bool, optional
        If ``True``, disable caching requests, defaults to False.

    Returns
    -------
    xarray.DataArray
        Elevations of the input coordinates as a ``xarray.DataArray``.
    """
    bbox = (min(xcoords), min(ycoords), max(xcoords), max(ycoords))

    ratio_min = 0.01
    ratio_x = abs((bbox[2] - bbox[0]) / bbox[0])
    ratio_y = abs((bbox[3] - bbox[1]) / bbox[1])
    if (ratio_x < ratio_min) or (ratio_y < ratio_min):
        rad = ratio_min * abs(bbox[0])
        bbox = (bbox[0] - rad, bbox[1] - rad, bbox[2] + rad, bbox[3] + rad)

    wms = WMS(
        ServiceURL().wms.nm_3dep,
        layers="3DEPElevation:None",
        outformat="image/tiff",
        crs=DEF_CRS,
        expire_after=expire_after,
        disable_caching=disable_caching,
        validation=False,
    )
    r_dict = wms.getmap_bybox(bbox, resolution, box_crs=crs)
    dem = utils.reproject_gtiff(r_dict, crs)

    if depression_filling:
        dem = utils.fill_depressions(dem)

    dem = dem.interp(x=list(xcoords), y=list(ycoords))
    valid_layers = wms.get_validlayers()
    return utils.rename_layers(dem, list(valid_layers))


class ElevationByCoords(BaseModel):
    """Elevation model class.

    Parameters
    ----------
    crs : str, optional
        Coordinate reference system of the input coordinates, defaults to ``EPSG:4326``.
    coords : list of tuple
        List of coordinates.
    source : str, optional
        Elevation source, defaults to ``tep``. Valid sources are: ``tnm``, ``tep``,
        and ``airmap``.
    expire_after : int, optional
        Expiration time for response caching in seconds, defaults to -1 (never expire).
    disable_caching : bool, optional
        If ``True``, disable caching requests, defaults to False.
    """

    crs: str = DEF_CRS
    coords: List[Tuple[float, float]]
    source: str = "tep"
    expire_after: float = EXPIRE
    disable_caching: bool = False

    @validator("crs")
    def _valid_crs(cls, v: str) -> str:
        return ogc_utils.validate_crs(v)

    @validator("coords")
    def _validate_coords(
        cls, v: List[Tuple[float, float]], values: Dict[str, str]
    ) -> gpd.GeoSeries:
        return gpd.GeoSeries(gpd.points_from_xy(*zip(*v)), crs=values["crs"])

    @validator("source")
    def _validate_source(cls, v: str) -> str:
        valid_sources = ["tnm", "airmap", "tep"]
        if v not in valid_sources:
            raise InvalidInputValue("source", valid_sources)
        return v

    def airmap(self) -> List[float]:
        """Return list of elevations in meters."""
        pts = self.coords.to_crs(DEF_CRS)  # type: ignore
        coords_chunks = tlz.partition_all(100, zip(pts.x, pts.y))
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
        elevations = list(
            tlz.pluck(
                "data",
                ar.retrieve_json(
                    urls, kwds, expire_after=self.expire_after, disable=self.disable_caching
                ),
            )
        )
        return list(tlz.concat(elevations))

    def tnm(self) -> List[float]:
        """Return list of elevations in meters."""
        pts = self.coords.to_crs(DEF_CRS)  # type: ignore
        kwds = [
            {"params": {"units": "Meters", "output": "json", "x": f"{x:.5f}", "y": f"{y:.5f}"}}
            for x, y in zip(pts.x, pts.y)
        ]
        resp = ar.retrieve_json(
            [ServiceURL().restful.nm_pqs] * len(kwds),
            kwds,
            max_workers=5,
            expire_after=self.expire_after,
            disable=self.disable_caching,
        )
        return [
            r["USGS_Elevation_Point_Query_Service"]["Elevation_Query"]["Elevation"] for r in resp
        ]

    def tep(self) -> List[float]:
        """Get elevation from 3DEP."""
        wms_3dep = WMS(
            ServiceURL().wms.nm_3dep,
            layers="3DEPElevation:None",
            outformat="image/tiff",
            crs="EPSG:3857",
            expire_after=self.expire_after,
            disable_caching=self.disable_caching,
            validation=False,
        )
        get_dem = tlz.partial(wms_3dep.getmap_bybox, resolution=10, box_crs=wms_3dep.crs)

        points_proj = self.coords.to_crs(wms_3dep.crs)  # type: ignore
        bounds = points_proj.buffer(30, cap_style=3)
        da_list = [geoutils.gtiff2xarray(get_dem(b.bounds)) for b in bounds]

        def get_value(da: xr.DataArray, x: float, y: float) -> float:
            nodata = da.attrs["nodatavals"][0]
            value = da.fillna(nodata).interp(x=[x], y=[y])
            return float(value.values[0, 0])

        return [get_value(da, p.x, p.y) for da, p in zip(da_list, points_proj)]


def elevation_bycoords(
    coords: List[Tuple[float, float]],
    crs: Union[str, pyproj.CRS] = DEF_CRS,
    source: str = "tep",
    expire_after: float = EXPIRE,
    disable_caching: bool = False,
) -> List[float]:
    """Get elevation for a list of coordinates.

    Parameters
    ----------
    coords : list of tuple
        Coordinates of target location as list of tuples ``[(x, y), ...]``.
    crs : str or pyproj.CRS, optional
        Spatial reference (CRS) of coords, defaults to ``EPSG:4326``.
    source : str, optional
        Data source to be used, default to ``airmap``. Supported sources are
        ``airmap`` (30 m resolution), ``tnm`` (using The National Map's Bulk Point
        Query Service with 10 m resolution) and ``tep`` (using 3DEP's WMS service
        at 10 m resolution). The ``tnm`` and ``tep`` sources are more accurate since they
        use the 1/3 arc-second DEM layer from 3DEP service but it is limited to the US.
        They both tend to be slower than the Airmap service. Note that ``tnm`` is bit unstable.
        It's recommended to use ``tep`` unless 10-m resolution accuracy is not necessary which
        in that case ``airmap`` is more appropriate.
    expire_after : int, optional
        Expiration time for response caching in seconds, defaults to -1 (never expire).
    disable_caching : bool, optional
        If ``True``, disable caching requests, defaults to False.

    Returns
    -------
    list of float
        Elevation in meter.
    """
    _crs = crs.to_string() if isinstance(crs, pyproj.CRS) else crs
    service = ElevationByCoords(
        crs=_crs,
        coords=coords,
        source=source,
        expire_after=expire_after,
        disable_caching=disable_caching,
    )
    return service.__getattribute__(source)()  # type: ignore


def check_3dep_availability(
    bbox: Tuple[float, float, float, float], crs: Union[str, pyproj.CRS] = DEF_CRS
) -> Dict[str, bool]:
    """Query 3DEP's resolution availability within a bounding box.

    This function checks availability of 3DEP's at the following resolutions:
    1 m, 3 m, 5 m, 10 m, 30 m, and 60 m.

    Parameters
    ----------
    bbox : tuple
        Bounding box as tuple of ``(min_x, min_y, max_x, max_y)``.
    crs : str or pyproj.CRS, optional
        Spatial reference (CRS) of bbox, defaults to ``EPSG:4326``.

    Returns
    -------
    dict
        True if bbox intersects 3DEP elevation for each available resolution.
        Keys are the supported resolutions and values are their availability.

    Examples
    --------
    >>> import py3dep
    >>> bbox = (-69.77, 45.07, -69.31, 45.45)
    >>> py3dep.check_3dep_availability(bbox)
    {'1m': True, '3m': False, '5m': False, '10m': True, '30m': True, '60m': False}
    """
    if not isinstance(bbox, Sequence) or len(bbox) != 4:
        raise InvalidInputType("bbox", "a tuple of 4 elements")

    res_layers = {
        "1m": 18,
        "3m": 19,
        "5m": 20,
        "10m": 21,
        "30m": 22,
        "60m": 23,
    }
    _bbox = ogc_utils.match_crs(bbox, crs, DEF_CRS)
    url = ServiceURL().restful.nm_3dep_index

    def _check(lyr: int) -> bool:
        wms = ArcGISRESTful(url, lyr)
        with contextlib.suppress(ServiceError, TypeError):
            return len(list(wms.oids_bygeom(_bbox))[0]) > 0
        return False

    return {res: _check(lyr) for res, lyr in res_layers.items()}
