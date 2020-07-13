"""Get data from 3DEP database."""
from pathlib import Path
from typing import List, Optional, Tuple, Union

import numpy as np
import xarray as xr
from shapely.geometry import Polygon

from . import utils
from .exceptions import InvalidInputType
from .utils import MatchCRS


def get_map(
    layers: Union[str, List[str]],
    geometry: Union[Polygon, Tuple[float, float, float, float]],
    resolution: float,
    geo_crs: str = "epsg:4326",
    crs: str = "epsg:4326",
    fill_holes: bool = False,
    data_dir: Optional[Union[str, Path]] = None,
) -> xr.DataArray:
    """Access to `3DEP <https://www.usgs.gov/core-science-systems/ngp/3dep>`__ service.

    The 3DEP service has multi-resolution sources so depending on the user
    provided resolution the data is resampled on server-side based
    on all the available data sources. The following layers are available:
    - "DEM"
    - "Hillshade Gray"
    - "Aspect Degrees"
    - "Aspect Map"
    - "GreyHillshade_elevationFill"
    - "Hillshade Multidirectional"
    - "Slope Map"
    - "Slope Degrees"
    - "Hillshade Elevation Tinted"
    - "Height Ellipsoidal"
    - "Contour 25"
    - "Contour Smoothed 25"

    Parameters
    ----------
    layer : str or list
        A valid 3DEP layer or a list of them
    geometry : shapely.geometry.Polygon
        A shapely Polygon in WGS 84 (epsg:4326).
    resolution : float
        The data resolution in meters. The width and height of the output are computed in pixel
        based on the geometry bounds and the given resolution.
    geo_crs : str, optional
        The spatial reference system of the input geometry, defaults to
        epsg:4326.
    crs : str, optional
        The spatial reference system to be used for requesting the data, defaults to
        epsg:4326.
    data_dir : str or Path, optional
        The directory to save the downloaded images, defaults to None which will only return
        the data as ``xarray.Dataset`` and doesn't save them as ``tiff`` images.
    fill_holes : bool, optional
        Whether to fill the holes in the geometry's interior, defaults to False.

    Returns
    -------
    xarray.DataArray
        The requeted data within the geometry
    """
    if not isinstance(geometry, (Polygon, tuple)):
        raise InvalidInputType("geometry", "Polygon or tuple of length 4")

    if isinstance(geometry, Polygon):
        _geometry = Polygon(geometry.exterior) if fill_holes else geometry
        _geometry = MatchCRS.geometry(_geometry, geo_crs, crs)
        bounds = _geometry.bounds
    else:
        utils.check_bbox(geometry)
        _geometry = MatchCRS.bounds(geometry, geo_crs, crs)
        bounds = _geometry

    r_dict = utils.wms_bybox(layers, bounds, resolution, box_crs=crs, crs=crs,)

    return utils.wms_toxarray(r_dict, _geometry, data_dir)


def elevation_bygrid(
    gridxy: List[List[float]], crs: str = "epsg:4326", resolution: float = 10
) -> np.ndarray:
    """Get elevation from DEM data for a list of coordinates.

    This function is intended for getting elevations for a gridded dataset.

    Parameters
    ----------
    gridxy : list of two lists of floats
        A list containing x- and y-coordinates of a mesh, [[x-coords], [y-coords]].
    crs : str, optional
        The spatial reference system of the input coords, defaults to epsg:4326.
    resolution : float
        The accuracy of the output, defaults to 10 m which is the highest
        available resolution that covers CONUS. Note that higher resolution
        increases computation time so chose this value with caution.

    Returns
    -------
    xarray.DataArray
        A data array with dims ``x`` and ``y``
    """
    gx, gy = gridxy
    bbox = (min(gx), min(gy), max(gx), max(gy))

    ratio_min = 0.01
    ratio_x = abs((bbox[2] - bbox[0]) / bbox[0])
    ratio_y = abs((bbox[3] - bbox[1]) / bbox[1])
    if (ratio_x < ratio_min) or (ratio_y < ratio_min):
        rad = ratio_min * abs(bbox[0])
        bbox = (bbox[0] - rad, bbox[1] - rad, bbox[2] + rad, bbox[3] + rad)

    dem = get_map("DEM", bbox, resolution, geo_crs=crs, crs=crs)

    elev = dem.interp(x=gx, y=gy)
    elev.name = "elevation"
    elev.attrs["units"] = "meters"
    return elev


def elevation_byloc(coord: Tuple[float, float], crs: str = "epsg:4326"):
    """Get elevation from USGS 3DEP service for a coordinate.

    Parameters
    ----------
    coord : tuple
        Coordinates of the location as a tuple
    crs : str, optional
        The spatial reference of the coord arg, defaults to EPSG:4326
    Returns
    -------
    float
        Elevation in meter
    """
    if not isinstance(coord, tuple) and len(tuple) != 2:
        raise InvalidInputType("coord", "tuple of length 2")

    lon, lat = MatchCRS.point(coord, crs, "epsg:4326")

    url = "https://nationalmap.gov/epqs/pqs.php"
    payload = {"output": "json", "x": lon, "y": lat, "units": "Meters"}
    r = utils.RetrySession().get(url, payload)
    root = r.json()["USGS_Elevation_Point_Query_Service"]
    elevation = float(root["Elevation_Query"]["Elevation"])
    if abs(elevation - (-1000000)) < 1e-3:
        raise ValueError(
            f"The elevation of the requested coordinate ({coord[0]}, {coord[1]}) cannot be found."
        )

    return elevation
