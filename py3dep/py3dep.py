"""Get data from 3DEP database."""
from typing import List, Tuple, Union

import numpy as np
import pygeoutils as geoutils
import rasterio as rio
import rasterio.warp as rio_warp
import xarray as xr
from pygeoogc import WMS, MatchCRS, RetrySession, ServiceURL
from shapely.geometry import Polygon

from .exceptions import InvalidInputType

DEF_CRS = "epsg:4326"


def get_map(
    layers: Union[str, List[str]],
    geometry: Union[Polygon, Tuple[float, float, float, float]],
    resolution: float,
    geo_crs: str = DEF_CRS,
    crs: str = DEF_CRS,
    fill_holes: bool = False,
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
    layers : str or list
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
    fill_holes : bool, optional
        Whether to fill the holes in the geometry's interior, defaults to False.

    Returns
    -------
    xarray.DataArray
        The requeted data within the geometry
    """
    if not isinstance(geometry, (Polygon, tuple)):
        raise InvalidInputType("geometry", "Polygon or tuple of length 4")

    _geometry = geoutils.geo2polygon(geometry, geo_crs, crs)
    _geometry = Polygon(_geometry.exterior) if fill_holes else _geometry

    _layers = layers if isinstance(layers, list) else [layers]
    if "DEM" in _layers:
        _layers[_layers.index("DEM")] = "None"

    _layers = [f"3DEPElevation:{lyr}" for lyr in _layers]

    wms = WMS(ServiceURL().wms.nm_3dep, layers=_layers, outformat="image/tiff", crs=crs)
    r_dict = wms.getmap_bybox(_geometry.bounds, resolution, box_crs=crs)

    return geoutils.gtiff2xarray(r_dict, _geometry, crs)


def elevation_bygrid(
    gridxy: Tuple[List[float], List[float]],
    crs: str,
    resolution: float,
    resampling: rio_warp.Resampling = rio_warp.Resampling.bilinear,
) -> xr.DataArray:
    """Get elevation from DEM data for a list of coordinates.

    This function is intended for getting elevations for a gridded dataset.

    Parameters
    ----------
    gridxy : tuple of two lists of floats
        A list containing x- and y-coordinates of a mesh, [[x-coords], [y-coords]].
    crs : str
        The spatial reference system of the input grid, defaults to epsg:4326.
    resolution : float
        The accuracy of the output, defaults to 10 m which is the highest
        available resolution that covers CONUS. Note that higher resolution
        increases computation time so chose this value with caution.
    resampling : rasterio.warp.Resampling
        The reasmpling method to use if the input crs is not in the supported
        3DEP's CRS list which are epsg:4326 and epsg:3857. It defaults to bilinear.
        The available methods can be found `here <https://rasterio.readthedocs.io/en/latest/api/rasterio.enums.html#rasterio.enums.Resampling>`__

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

    req_crs = crs if crs.lower() in [DEF_CRS, "epsg:3857"] else DEF_CRS

    wms = WMS(
        ServiceURL().wms.nm_3dep, layers="3DEPElevation:None", outformat="image/tiff", crs=req_crs,
    )

    r_dict = wms.getmap_bybox(bbox, resolution, box_crs=crs,)

    def reproject(content):
        with rio.MemoryFile() as memfile:
            memfile.write(content)
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
                                resampling=resampling,
                            )

                    da = xr.open_rasterio(vrt)
                    try:
                        da = da.squeeze("band", drop=True)
                    except ValueError:
                        pass

                    da.name = "elevation"
                    da.attrs["transform"] = transform
                    da.attrs["res"] = (transform[0], transform[4])
                    da.attrs["bounds"] = tuple(vrt.bounds)
                    da.attrs["nodatavals"] = vrt.nodatavals
                    da.attrs["crs"] = vrt.crs
        return da

    _elev = xr.merge([reproject(c) for c in r_dict.values()])
    elev = _elev.elevation.interp(x=gx, y=gy)
    elev.attrs["units"] = "meters"
    return elev


def elevation_byloc(coord: Tuple[float, float], crs: str = DEF_CRS):
    """Get elevation from USGS 3DEP service for a coordinate.

    Parameters
    ----------
    coord : tuple
        Coordinates of the location as a tuple
    crs : str, optional
        The spatial reference of the input coord, defaults to epsg:4326 (lon, lat)

    Returns
    -------
    float
        Elevation in meter
    """
    if not isinstance(coord, tuple) or len(coord) != 2:
        raise InvalidInputType("coord", "tuple of length 2", "(x, y)")

    lon, lat = MatchCRS.coords(([coord[0]], [coord[1]]), crs, DEF_CRS)

    url = "https://nationalmap.gov/epqs/pqs.php"
    payload = {"output": "json", "x": lon[0], "y": lat[0], "units": "Meters"}
    r = RetrySession().get(url, payload)
    root = r.json()["USGS_Elevation_Point_Query_Service"]
    elevation = float(root["Elevation_Query"]["Elevation"])

    if abs(elevation - (-1000000)) < 1e-3:
        raise ValueError(
            f"The elevation of the requested coordinate ({coord[0]}, {coord[1]}) cannot be found."
        )

    return elevation


def deg2mpm(da: xr.DataArray) -> xr.DataArray:
    """Convert ``xarray.Data[Array,set]`` from degree to meter/meter."""
    attrs = da.attrs
    da = np.tan(np.deg2rad(da))
    da.attrs = attrs
    da.name = "slope"
    da.attrs["units"] = "meters/meters"
    return da
