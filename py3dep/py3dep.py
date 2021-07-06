"""Get data from 3DEP database."""
from itertools import product
from typing import Iterator, List, Optional, Tuple, Union

import async_retriever as ar
import cytoolz as tlz
import numpy as np
import pygeoutils as geoutils
import rasterio as rio
import rasterio.warp as rio_warp
import xarray as xr
from affine import Affine
from pygeoogc import WMS, ServiceURL, utils
from shapely.geometry import MultiPolygon, Polygon

try:
    import richdem as rd
except ImportError:
    rd = None
from .exceptions import InvalidInputType, InvalidInputValue

DEF_CRS = "epsg:4326"


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
        ``EPSG:4326``.

    Returns
    -------
    dict
        A dict where the keys are the layer name and values are the returned response
        from the WMS service as bytes. You can use ``utils.create_dataset`` function
        to convert the responses to ``xarray.Dataset``.
    """
    _geometry = geoutils.geo2polygon(geometry, geo_crs, crs)

    _layers = layers if isinstance(layers, list) else [layers]
    if "DEM" in _layers:
        _layers[_layers.index("DEM")] = "None"

    _layers = [f"3DEPElevation:{lyr}" for lyr in _layers]

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
    dim_names: Optional[Tuple[str, str]] = None,
    depression_filling: bool = False,
    resampling: rio_warp.Resampling = rio_warp.Resampling.bilinear,
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
    dim_names : tuple
        A tuple of length two containing the coordinate names, defaults to ``("x", "y")``.
    depression_filling : bool, optional
        Fill depressions before sampling using
        `RichDEM <https://richdem.readthedocs.io/en/latest/>`__ package, defaults to False.
    resampling : rasterio.warp.Resampling
        The reasmpling method to use if the input crs is not in the supported
        3DEP's CRS list which are ``EPSG:4326`` and ``EPSG:3857``.
        It defaults to bilinear. The available methods can be found
        `here <https://rasterio.readthedocs.io/en/latest/api/rasterio.enums.html#rasterio.enums.Resampling>`__

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
    elev_arr = _sample_tiff(
        r_dict["3DEPElevation:None_dd_0_0"], coords, crs, depression_filling, resampling
    )

    return xr.DataArray(
        elev_arr.reshape((len(xcoords), len(ycoords))),
        dims=dim_names,
        coords=[xcoords, ycoords],
        name="elevation",
        attrs={"units": "meters"},
    )


def fill_depressions(
    dem: Union[xr.DataArray, np.ndarray],
    nodata: Optional[float] = None,
    crs: Optional[str] = None,
    transform: Optional[Affine] = None,
) -> xr.DataArray:
    """Fill depressions and adjust flat areas in a DEM using RichDEM.

    Parameters
    ----------
    dem : xarray.DataArray or numpy.ndarray
        Digital Elevation Model.
    nodata : float, optional
        Value for cell with no data, defaults to None. If input DEM is an ``xarray.DataArray``
        nodata will be set to ``dem.nodatavals[0]``.
    crs : str, optional
        Coordinate reference system for the input DEM, defaults to None. If input DEM
         is an ``xarray.DataArray`` crs will be set to ``dem.crs``.
    transform : str, optional
        Transform of the input DEM, defaults to None. If input DEM is an ``xarray.DataArray``
        transform will be set to ``dem.transform``.

    Returns
    -------
    xarray.DataArray
        Conditioned DEM.
    """
    if rd is None:
        raise ImportError(
            "Depression filling operation uses richdem package which is not installed."
        )

    attrs = {"nodata": nodata, "crs": crs, "transform": transform}
    if isinstance(dem, xr.DataArray):
        try:
            attrs = {"nodata": dem.nodatavals[0], "crs": dem.crs, "transform": dem.transform}
        except AttributeError:
            raise AttributeError("Input DEM must have nodatavals, crs, and transform attributes.")

    rda = rd.rdarray(dem, no_data=attrs["nodata"])
    rda.projection = attrs["crs"]
    rda.geotransform = attrs["transform"]
    rda = rd.FillDepressions(rda, epsilon=False)
    rda = rd.ResolveFlats(rda)
    if isinstance(dem, xr.DataArray):
        dem.data = rda
        return dem
    return rda


def _sample_tiff(
    content: bytes,
    coords: Union[List[Tuple[float, float]], Iterator[Tuple[float, float]]],
    crs: str,
    depression_filling: bool,
    resampling: rio_warp.Resampling,
) -> np.ndarray:
    """Sample a tiff response for a list of coordinates.

    Parameters
    ----------
    content : bytes
    coords : list of tuples
        A list containing x- and y-coordinates of a mesh, ``[(x, y), ...]``.
    crs : str
        The spatial reference system of the input grid, defaults to epsg:4326.
    depression_filling : bool, optional
        Fill depressions before sampling using
        `RichDEM <https://richdem.readthedocs.io/en/latest/>`__ package, defaults to False.
    resampling : rasterio.warp.Resampling
        The reasmpling method to use if the input crs is not in the supported
        3DEP's CRS list which are epsg:4326 and ``epsg:3857``. It defaults to bilinear.
        The available methods can be found
        `here <https://rasterio.readthedocs.io/en/latest/api/rasterio.enums.html#rasterio.enums.Resampling>`__

    Returns
    -------
    numpy.ndarray
        An array of elevations where its index matches the input coords list.
    """
    with rio.MemoryFile() as memfile:
        memfile.write(content)
        with memfile.open() as src:
            transform, width, height = rio_warp.calculate_default_transform(
                src.crs, crs, src.width, src.height, *src.bounds
            )
            kwargs = src.meta.copy()
            kwargs.update({"crs": crs, "transform": transform, "width": width, "height": height})
            if depression_filling:
                src_nodata, src_crs = geoutils.get_nodata_crs(content)
                data = fill_depressions(src.read(1), src_nodata, src_crs, src.transform)
            else:
                data = rio.band(src, 1)
            with rio.vrt.WarpedVRT(src, **kwargs) as vrt:
                if crs != src.crs:
                    rio_warp.reproject(
                        source=data,
                        destination=rio.band(vrt, 1),
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
    """Get elevation within a bounding box.

    This function is intended for getting elevations for a gridded dataset.

    Parameters
    ----------
    bbox : tuple of two lists of floats
        Bounding box as a tuple of length of 4 ``(west, south, east, north)``.
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


def elevation_bycoords(
    coords: List[Tuple[float, float]], crs: str = DEF_CRS, source: str = "airmap"
) -> List[float]:
    """Get elevation from Airmap at 1-arc resolution (~30 m) for a list of coordinates.

    Parameters
    ----------
    coords : list of tuples
        Coordinates of target location as list of tuples ``[(x, y), ...]``.
    crs : str, optional
        Spatial reference (CRS) of coords, defaults to ``EPSG:4326``.
    source : str, optional
        Data source to be used, default to ``airmap``. Supported sources are
        ``airmap`` (30 m resolution) and ``tnm`` (using The National Map's Bulk Point
        Query Service). The ``tnm`` source is more accurate since it uses the highest
        available resolution DEM automatically but it is limited to the US and tends to
        be slower than ``airmap``.

    Returns
    -------
    list of float
        Elevation in meter.
    """
    coords = utils.match_crs(coords, crs, DEF_CRS)

    valid_sources = ["airmap", "tnm"]
    if source not in valid_sources:
        raise InvalidInputValue("source", valid_sources)

    if source == "airmap":
        coords_chunks = tlz.partition_all(100, coords)
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
            for lon, lat in coords
        )
    )
    resp = ar.retrieve(urls, "json", kwds, max_workers=500)
    return [r["USGS_Elevation_Point_Query_Service"]["Elevation_Query"]["Elevation"] for r in resp]


def deg2mpm(da: xr.DataArray) -> xr.DataArray:
    """Convert slope from degree to meter/meter."""
    attrs = da.attrs
    da = np.tan(np.deg2rad(da))
    da.attrs = attrs
    da.name = "slope"
    da.attrs["units"] = "meters/meters"
    return da
