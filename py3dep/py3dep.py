"""Get data from 3DEP database."""
from typing import List, Tuple, Union

import async_retriever as ar
import cytoolz as tlz
import pygeoutils as geoutils
import pyproj
import xarray as xr
from pydantic import BaseModel, validator
from pygeoogc import WMS, InvalidInputType, InvalidInputValue, ServiceURL
from pygeoogc import utils as ogc_utils
from shapely.geometry import MultiPolygon, Polygon

from . import utils

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
__all__ = ["get_map", "elevation_bygrid", "elevation_bycoords"]


def get_map(
    layers: Union[str, List[str]],
    geometry: Union[Polygon, MultiPolygon, Tuple[float, float, float, float]],
    resolution: float,
    geo_crs: str = DEF_CRS,
    crs: str = DEF_CRS,
) -> Union[xr.DataArray, xr.Dataset]:
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

    Returns
    -------
    xarray.DataArray or xarray.Dataset
        The requested topographic data as an xarray.DataArray or xarray.Dataset.
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
        ServiceURL().wms.nm_3dep, layers="3DEPElevation:None", outformat="image/tiff", crs=DEF_CRS
    )
    r_dict = wms.getmap_bybox(bbox, resolution, box_crs=crs)
    dem = utils.reproject_gtiff(r_dict, crs)

    if depression_filling:
        dem = utils.fill_depressions(dem)

    return dem.interp(x=xcoords, y=ycoords)


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
    source: str = "airmap"

    @validator("crs")
    def _valid_crs(cls, v):
        try:
            return pyproj.CRS(v)
        except pyproj.exceptions.CRSError as ex:
            raise InvalidInputType("crs", "a valid CRS") from ex

    @validator("coords")
    def _validate_coords(cls, v, values):
        return ogc_utils.match_crs(v, values["crs"], DEF_CRS)

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
        headers = {
            "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8",
            "Accept-Language": "en-US,en;q=0.5",
            "Accept-Encoding": "gzip,deflate,br",
            "Connection": "keep-alive",
            "Upgrade-Insecure-Requests": "1",
            "Sec-Fetch-Dest": "document",
            "Sec-Fetch-Mode": "navigate",
            "Sec-Fetch-Site": "same-origin",
            "Sec-Fetch-User": "?1",
        }
        urls, kwds = zip(
            *(
                (
                    ServiceURL().restful.nm_pqs,
                    {
                        "params": {
                            "units": "Meters",
                            "output": "json",
                            "y": f"{lat:.5f}",
                            "x": f"{lon:.5f}",
                        },
                        "headers": headers,
                    },
                )
                for lon, lat in self.coords
            )
        )
        resp = ar.retrieve(urls, "json", kwds, max_workers=5)
        return [
            r["USGS_Elevation_Point_Query_Service"]["Elevation_Query"]["Elevation"] for r in resp
        ]


def elevation_bycoords(
    coords: List[Tuple[float, float]], crs: Union[str, pyproj.CRS] = DEF_CRS, source: str = "airmap"
) -> List[float]:
    """Get elevation from Airmap at 1-arc resolution (~30 m) for a list of coordinates.

    Parameters
    ----------
    coords : list of tuples
        Coordinates of target location as list of tuples ``[(x, y), ...]``.
    crs : str or pyproj.CRS, optional
        Spatial reference (CRS) of coords, defaults to ``EPSG:4326``.
    source : str, optional
        Data source to be used, default to ``airmap``. Supported sources are
        ``airmap`` (30 m resolution) and ``tnm`` (using The National Map's Bulk Point
        Query Service with 10 m resolution). The ``tnm`` source is more accurate since it
        uses the 1/3 arc-second DEM layer from 3DEP service but it is limited to the US.
        It also tends to be slower than the Airmap service and more unstable.
        It's recommended to use ``airmap`` unless you need 10-m resolution accuracy.

    Returns
    -------
    list of float
        Elevation in meter.
    """
    _crs = crs.to_string() if isinstance(crs, pyproj.CRS) else crs
    service = ElevationByCoords(crs=_crs, coords=coords, source=source)
    return getattr(service, source)()
