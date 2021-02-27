"""Get data from 3DEP database."""
from itertools import product
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Tuple, Union

import cytoolz as tlz
import numpy as np
import pygeoutils as geoutils
import rasterio as rio
import rasterio.warp as rio_warp
import xarray as xr
from pygeoogc import WMS, MatchCRS, RetrySession, ServiceURL
from shapely.geometry import MultiPolygon, Polygon

from .exceptions import InvalidInputType

DEF_CRS = "epsg:4326"


def get_map(
    layers: Union[str, List[str]],
    geometry: Union[Polygon, Tuple[float, float, float, float]],
    resolution: float,
    geo_crs: str = DEF_CRS,
    crs: str = DEF_CRS,
    output_dir: Optional[Union[str, Path]] = None,
) -> Dict[str, bytes]:
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
    geometry : Polygon, MultiPolygon, or tuple
        A shapely Polygon or a bounding box (west, south, east, north)
    resolution : float
        The data resolution in meters. The width and height of the output are computed in pixel
        based on the geometry bounds and the given resolution.
    geo_crs : str, optional
        The spatial reference system of the input geometry, defaults to
        epsg:4326.
    crs : str, optional
        The spatial reference system to be used for requesting the data, defaults to
        epsg:4326.
    output_dir : str or Path, optional
        The output directory to also save the map as GTiff file(s), defaults to None.

    Returns
    -------
    dict
        A dict where the keys are the layer name and values are the returned response
        from the WMS service as bytes. You can use ``utils.create_dataset`` function
        to convert the responses to ``xarray.Dataset``.
    """
    if not isinstance(geometry, (Polygon, MultiPolygon, tuple)):
        raise InvalidInputType("geometry", "Polygon or tuple of length 4")

    _geometry = geoutils.geo2polygon(geometry, geo_crs, crs)

    _layers = layers if isinstance(layers, list) else [layers]
    if "DEM" in _layers:
        _layers[_layers.index("DEM")] = "None"

    _layers = [f"3DEPElevation:{lyr}" for lyr in _layers]

    wms = WMS(ServiceURL().wms.nm_3dep, layers=_layers, outformat="image/tiff", crs=crs)
    r_dict = wms.getmap_bybox(_geometry.bounds, resolution, box_crs=crs)

    if output_dir:
        geoutils.gtiff2file(r_dict, _geometry, crs, output_dir)

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
    dim_names: Optional[Tuple[str, str]] = None,
    resampling: rio_warp.Resampling = rio_warp.Resampling.bilinear,
) -> xr.DataArray:
    """Get elevation from DEM data for a grid.

    This function is intended for getting elevations for a gridded dataset.

    Parameters
    ----------
    xcoords : tuple of two lists of floats
        A list containing x-coordinates of a mesh.
    ycoords : tuple of two lists of floats
        A list containing y-coordinates of a mesh.
    crs : str
        The spatial reference system of the input grid, defaults to epsg:4326.
    resolution : float
        The accuracy of the output, defaults to 10 m which is the highest
        available resolution that covers CONUS. Note that higher resolution
        increases computation time so chose this value with caution.
    dim_names : tuple
        A tuple of length two containing the coordinate names, defaults to ["x", "y"]
    resampling : rasterio.warp.Resampling
        The reasmpling method to use if the input crs is not in the supported
        3DEP's CRS list which are epsg:4326 and epsg:3857. It defaults to bilinear.
        The available methods can be found `here <https://rasterio.readthedocs.io/en/latest/api/rasterio.enums.html#rasterio.enums.Resampling>`__

    Returns
    -------
    xarray.DataArray
        An data array with name elevation and the given dim names.
    """
    if dim_names is None:
        dim_names = ("x", "y")

    bbox = (min(xcoords), min(ycoords), max(xcoords), max(ycoords))
    r_dict = _elevation_bybox(bbox, crs, resolution)
    coords = product(xcoords, ycoords)
    elev_arr = _sample_tiff(r_dict["3DEPElevation:None_dd_0_0"], coords, crs, resampling)

    return xr.DataArray(
        elev_arr.reshape((len(xcoords), len(ycoords))),
        dims=dim_names,
        coords=[xcoords, ycoords],
        name="elevation",
        attrs={"units": "meters"},
    )


def _sample_tiff(
    content: bytes,
    coords: Union[List[Tuple[float, float]], Iterator[Tuple[float, float]]],
    crs: str,
    resampling: rio_warp.Resampling,
) -> np.ndarray:
    """Sample a tiff response for a list of coordinates.

    Parameters
    ----------
    content : bytes
    coords : list of tuples
        A list containing x- and y-coordinates of a mesh, [(x, y), ...].
    crs : str
        The spatial reference system of the input grid, defaults to epsg:4326.
    resolution : float
    resampling : rasterio.warp.Resampling
        The reasmpling method to use if the input crs is not in the supported
        3DEP's CRS list which are epsg:4326 and epsg:3857. It defaults to bilinear.
        The available methods can be found `here <https://rasterio.readthedocs.io/en/latest/api/rasterio.enums.html#rasterio.enums.Resampling>`__

    Returns
    -------
    numpy.ndarray
        An array of elevations where its index matches the input coords list
    """
    with rio.MemoryFile() as memfile:
        memfile.write(content)
        with memfile.open() as src:
            transform, width, height = rio_warp.calculate_default_transform(
                src.crs, crs, src.width, src.height, *src.bounds
            )
            kwargs = src.meta.copy()
            kwargs.update({"crs": crs, "transform": transform, "width": width, "height": height})

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
                return np.array([e.item() for e in vrt.sample(coords)])


def _elevation_bybox(
    bbox: Tuple[float, float, float, float],
    crs: str,
    resolution: float,
) -> xr.DataArray:
    """Get elevation from DEM data for a list of coordinates.

    This function is intended for getting elevations for a gridded dataset.

    Parameters
    ----------
    bbox : tuple of two lists of floats
        A list containing x- and y-coordinates of a mesh, [[x-coords], [y-coords]].
    crs : str
        The spatial reference system of the input grid, defaults to epsg:4326.
    resolution : float
        The accuracy of the output, defaults to 10 m which is the highest
        available resolution that covers CONUS. Note that higher resolution
        increases computation time so chose this value with caution.

    Returns
    -------
    numpy.ndarray
        An array of elevations where its index matches the input gridxy list
    """
    if not isinstance(bbox, tuple) or len(bbox) != 4:
        raise InvalidInputType("bbox", "tuple of length 4")

    ratio_min = 0.01
    ratio_x = abs((bbox[2] - bbox[0]) / bbox[0])
    ratio_y = abs((bbox[3] - bbox[1]) / bbox[1])
    if (ratio_x < ratio_min) or (ratio_y < ratio_min):
        rad = ratio_min * abs(bbox[0])
        bbox = (bbox[0] - rad, bbox[1] - rad, bbox[2] + rad, bbox[3] + rad)

    req_crs = crs if crs.lower() in [DEF_CRS, "epsg:3857"] else DEF_CRS
    wms = WMS(
        ServiceURL().wms.nm_3dep, layers="3DEPElevation:None", outformat="image/tiff", crs=req_crs
    )
    return wms.getmap_bybox(bbox, resolution, box_crs=crs)


def elevation_bycoords(coords: List[Tuple[float, float]], crs: str = DEF_CRS) -> List[int]:
    """Get elevation from Airmap for a list of coordinates.

    Parameters
    ----------
    coords : list of tuples
        Coordinates of the location as a tuple
    crs : str, optional
        The spatial reference of the input coord, defaults to epsg:4326 (lon, lat)

    Returns
    -------
    list of int
        Elevation in meter
    """
    if not isinstance(coords, (list, Iterator)):
        raise InvalidInputType("coord", "list (or iterator) of tuples of length 2", "[(x, y), ...]")

    if isinstance(coords, list) and any(len(c) != 2 for c in coords):
        raise InvalidInputType("coord", "list of tuples of length 2", "[(x, y), ...]")

    coords_reproj = zip(*MatchCRS.coords(tuple(zip(*coords)), crs, DEF_CRS))
    coords_reproj = tlz.partition_all(100, coords_reproj)

    headers = {"Content-Type": "application/json", "charset": "utf-8"}
    elevations = []
    for chunk in coords_reproj:
        payload = {"points": ",".join(f"{lat},{lon}" for lon, lat in chunk)}
        resp = RetrySession().get(ServiceURL().restful.airmap, payload=payload, headers=headers)
        elevations.append(resp.json()["data"])

    return list(tlz.concat(elevations))


def deg2mpm(da: xr.DataArray) -> xr.DataArray:
    """Convert ``xarray.Data[Array,set]`` from degree to meter/meter."""
    attrs = da.attrs
    da = np.tan(np.deg2rad(da))
    da.attrs = attrs
    da.name = "slope"
    da.attrs["units"] = "meters/meters"
    return da
