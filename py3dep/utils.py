"""Some utilities for Py3DEP."""
import os
from concurrent import futures
from itertools import product
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Mapping, Optional, Tuple, Union, ValuesView

import numpy as np
import pyproj
import rasterio as rio
import xarray as xr
from owslib.wms import WebMapService
from rasterio import mask as rio_mask
from rasterio import warp as rio_warp
from requests import Response, Session
from requests.adapters import HTTPAdapter
from requests.exceptions import HTTPError, RequestException, RetryError, Timeout
from shapely.geometry import LineString, Point, Polygon, box, mapping, shape
from urllib3 import Retry

from .exceptions import InvalidInputType, InvalidInputValue


def threading(
    func: Callable,
    iter_list: Iterable,
    param_list: Optional[List[Any]] = None,
    max_workers: int = 8,
) -> List[Any]:
    """Run a function in parallel with threading.

    Notes
    -----
    This function is suitable for IO intensive functions.

    Parameters
    ----------
    func : function
        The function to be ran in threads
    iter_list : list
        The iterable for the function
    param_list : list, optional
        List of other parameters, defaults to an empty list
    max_workers : int, optional
        Maximum number of threads, defaults to 8

    Returns
    -------
    list
        A list of function returns for each iterable. The list is not ordered.
    """
    data = []
    param_list = [] if param_list is None else param_list
    with futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_itr = {executor.submit(func, itr, *param_list): itr for itr in iter_list}
        for future in futures.as_completed(future_to_itr):
            itr = future_to_itr[future]
            try:
                data.append(future.result())
            except Exception as exc:
                raise Exception(f"{itr}: {exc}")
    return data


def create_dataset(
    content: bytes,
    geometry: Union[Polygon, Tuple[float, float, float, float]],
    name: str,
    fpath: Optional[Union[str, Path]] = None,
) -> Union[xr.Dataset, xr.DataArray]:
    """Create dataset from a response clipped by a geometry.

    Parameters
    ----------
    content : requests.Response
        The response to be processed
    geometry : Polygon
        The geometry for masking the data
    name : str
        Variable name in the dataset
    fpath : str or Path, optinal
        The path save the file, defaults to None i.e., don't save as an image.

    Returns
    -------
    xarray.Dataset
        Generated xarray DataSet or DataArray
    """
    geom = [box(*geometry)] if isinstance(geometry, tuple) else [geometry]

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

            masked, transform = rio_mask.mask(src, geom, crop=True, nodata=nodata)
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


class MatchCRS:
    """Match CRS of a input geometry (Polygon, bbox, coord) with the output CRS.

    Parameters
    ----------
    geometry : tuple or Polygon
        The input geometry (Polygon, bbox, coord)
    in_crs : str
        The spatial reference of the input geometry
    out_crs : str
        The target spatial reference
    """

    @staticmethod
    def geometry(geom: Polygon, in_crs: str, out_crs: str) -> Polygon:
        if not isinstance(geom, Polygon):
            raise InvalidInputType("geom", "Polygon")

        return shape(rio_warp.transform_geom(in_crs, out_crs, mapping(geom)))

    @staticmethod
    def bounds(
        geom: Tuple[float, float, float, float], in_crs: str, out_crs: str
    ) -> Tuple[float, float, float, float]:
        if not isinstance(geom, tuple) and len(geom) != 4:
            raise InvalidInputType("geom", "tuple of length 4")

        return rio_warp.transform_bounds(in_crs, out_crs, *geom)

    @staticmethod
    def point(geom: Tuple[float, float], in_crs: str, out_crs: str) -> Tuple[float, float]:
        if not isinstance(geom, tuple) and len(geom) != 2:
            raise InvalidInputType("geom", "tuple of length 2")

        pts = rio_warp.transform(in_crs, out_crs, [geom[0]], [geom[1]])
        return (pts[0][0], pts[1][0])


def check_bbox(bbox: Tuple[float, float, float, float]) -> None:
    """Check if an input inbox is a tuple of length 4."""
    if not isinstance(bbox, tuple) or len(bbox) != 4:
        raise InvalidInputType("bbox", "tuple", "(west, south, east, north)")


def bbox_resolution(
    bbox: Tuple[float, float, float, float], resolution: float, bbox_crs: str = "epsg:4326"
) -> Tuple[int, int]:
    """Image size of a bounding box WGS84 for a given resolution in meters.

    Parameters
    ----------
    bbox : tuple
        A bounding box in WGS84 (west, south, east, north)
    resolution : float
        The resolution in meters
    bbox_crs : str, optional
        The spatial reference of the input bbox, default to EPSG:4326.

    Returns
    -------
    tuple
        The width and height of the image
    """
    check_bbox(bbox)

    bbox = MatchCRS.bounds(bbox, bbox_crs, "epsg:4326")
    west, south, east, north = bbox
    geod = pyproj.Geod(ellps="WGS84")

    linex = LineString([Point(west, south), Point(east, south)])
    delx = geod.geometry_length(linex)

    liney = LineString([Point(west, south), Point(west, north)])
    dely = geod.geometry_length(liney)

    return int(delx / resolution), int(dely / resolution)


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


def wms_bybox(
    layers: Union[str, List[str]],
    bbox: Tuple[float, float, float, float],
    resolution: float,
    box_crs: str = "epsg:4326",
    crs: str = "epsg:4326",
    version: str = "1.3.0",
) -> Dict[str, bytes]:
    """Get data from a 3DEP service within a geometry or bounding box.

    Parameters
    ----------
    layers : str or list
        A layer or a list of layers from the service to be downloaded. You can pass an empty
        string to get a list of available layers.
    box : tuple
        A bounding box for getting the data.
    resolution : float
        The output resolution in meters. The width and height of output are computed in pixel
        based on the geometry bounds and the given resolution.
    box_crs : str, optional
        The spatial reference system of the input bbox, defaults to
        epsg:4326.
    crs : str, optional
        The spatial reference system to be used for requesting the data, defaults to
        epsg:4326.
    version : str, optional
        The WMS service version which should be either 1.1.1 or 1.3.0, defaults to 1.3.0.

    Returns
    -------
    dict
        A dict where the keys are the layer name and values are the returned response
        from the WMS service as bytes.
    """
    url = "https://elevation.nationalmap.gov/arcgis/services/3DEPElevation/ImageServer/WMSServer"
    wms = WebMapService(url, version=version)

    valid_layers = [wms[lyr].name.split(":")[-1] for lyr in list(wms.contents)]
    del valid_layers[valid_layers.index("0")]
    valid_layers[valid_layers.index("None")] = "DEM"

    if not isinstance(layers, (str, list)):
        raise InvalidInputType("layers", "str or list")

    _layers = [layers] if isinstance(layers, str) else layers

    if any(lyr not in valid_layers for lyr in _layers):
        raise InvalidInputValue("layer", (lyr for lyr in valid_layers))

    _layers = ["None" if lyr == "DEM" else lyr for lyr in _layers]
    _layers = [f"3DEPElevation:{lyr}" for lyr in _layers]

    valid_crss = {lyr: [s.lower() for s in wms[lyr].crsOptions] for lyr in _layers}
    if any(crs not in valid_crss[lyr] for lyr in _layers):
        _valid_crss = (f"{lyr}: {', '.join(cs)}\n" for lyr, cs in valid_crss.items())
        raise InvalidInputValue("CRS", _valid_crss)

    check_bbox(bbox)
    _bbox = MatchCRS.bounds(bbox, box_crs, crs)
    width, height = bbox_resolution(_bbox, resolution, crs)
    bounds, widths = vsplit_bbox(_bbox, resolution, crs)
    _bounds = [(*bw[0], i, bw[1]) for i, bw in enumerate(zip(bounds, widths))]

    def getmap(args):
        lyr, bnds = args
        _bbox, res_count, _width = bnds[:-2], bnds[-2], bnds[-1]
        img = wms.getmap(
            layers=[lyr], srs=crs, bbox=_bbox, size=(_width, height), format="image/tiff",
        )
        return (f"{lyr}_{res_count}", img.read())

    return dict(threading(getmap, product(_layers, _bounds), max_workers=len(_layers)))


def vsplit_bbox(
    bbox: Tuple[float, float, float, float],
    resolution: float,
    box_crs: str = "epsg:4326",
    max_pixel: int = 8000000,
) -> Tuple[List[Tuple[float, float, float, float]], List[int]]:
    """Split the bounding box vertically for WMS requests.

    Parameters
    ----------
    bbox : tuple
        A bounding box; (west, south, east, north)
    resolution : float
        The target resolution for a WMS request in meters.
    box_crs : str, optional
        The spatial reference of the input bbox, default to EPSG:4326.
    max_pixel : int, opitonal
        The maximum allowable number of pixels (width x height) for a WMS requests,
        defaults to 8 million based on some trial-and-error.

    Returns
    -------
    tuple
        The first element is a list of bboxes and the second one is width of the last bbox
    """
    _crs = "epsg:4326"
    check_bbox(bbox)
    _bbox = MatchCRS.bounds(bbox, box_crs, _crs)
    width, height = bbox_resolution(_bbox, resolution, _crs)

    # Max number of cells that 3DEP can handle safely is about 8 mil.
    # We need to divide the domain into boxes with cell count of 8 mil.
    # We fix the height and incremenet the width.
    if (width * height) > max_pixel:
        _width = int(max_pixel / height)
        stepx = _width * resolution

        geod = pyproj.Geod(ellps="WGS84")

        west, south, east, north = _bbox
        az = geod.inv(west, south, east, south)[0]
        lons = [west]

        while lons[-1] < east:
            lons.append(geod.fwd(lons[-1], south, az, stepx)[0])

        lons[-1] = east
        bboxs = [(left, south, right, north) for left, right in zip(lons[:-1], lons[1:])]
        bboxs = [MatchCRS.bounds(i, _crs, box_crs) for i in bboxs]
        widths = [_width for _ in range(len(bboxs))]
        widths[-1] = width - (len(bboxs) - 1) * _width
    else:
        bboxs = [bbox]
        widths = [width]
    return bboxs, widths


def wms_toxarray(
    r_dict: Dict[str, bytes],
    geometry: Union[Polygon, Tuple[float, float, float, float]],
    data_dir: Optional[Union[str, Path]] = None,
) -> Union[xr.DataArray, xr.Dataset]:
    """Convert responses from ``wms_bybox`` to ``xarray.Dataset`` or ``xarray.DataAraay``.

    Parameters
    ----------
    r_dict : dict
        The output of ``wms_bybox`` function.
    geometry : Polygon or tuple
        The geometry to mask the data that should be in the same CRS as the r_dict.
    data_dir : str or Path, optional
        The directory to save the output as ``tiff`` images.

    Returns
    -------
    xarray.Dataset or xarray.DataAraay
        The dataset or data array based on the number of variables.
    """
    var_name = {
        lyr: lyr.split("_")[0].split(":")[-1].lower().replace(" ", "_") for lyr in r_dict.keys()
    }
    var_name = {lyr: "dem" if var == "none" else var for lyr, var in var_name.items()}

    if data_dir is not None:
        check_dir(Path(data_dir, "dummy"))

    _fpath = {lyr: f'{lyr.split(":")[-1].lower().replace(" ", "_")}.tiff' for lyr in r_dict.keys()}
    fpath = {lyr: None if data_dir is None else Path(data_dir, fn) for lyr, fn in _fpath.items()}

    ds = xr.merge(
        [create_dataset(r, geometry, var_name[lyr], fpath[lyr]) for lyr, r in r_dict.items()]
    )

    def mask_da(da):
        if not np.isnan(da.nodatavals[0]):
            _da = da.where(da < da.nodatavals[0], drop=True)
            _da.attrs["nodatavals"] = (np.nan,)
            return _da
        return da

    ds = ds.apply(mask_da)

    if len(ds.variables) - len(ds.dims) == 1:
        ds = ds[list(ds.keys())[0]]
    return ds


def deg2mpm(da: xr.DataArray) -> xr.DataArray:
    """Convert ``xarray.Data[Array,set]`` from degree to meter/meter."""
    attrs = da.attrs
    da = np.tan(np.deg2rad(da))
    da.attrs = attrs
    da.attrs["units"] = "meters/meters"
    return da


class RetrySession:
    """Configures the passed-in session to retry on failed requests.

    The fails can be due to connection errors, specific HTTP response
    codes and 30X redirections. The code is based on:
    https://github.com/bustawin/retry-requests

    Parameters
    ----------
    retries : int, optional
        The number of maximum retries before raising an exception, defaults to 5.
    backoff_factor : float, optional
        A factor used to compute the waiting time between retries, defaults to 0.5.
    status_to_retry : tuple, optional
        A tuple of status codes that trigger the reply behaviour, defaults to (500, 502, 504).
    prefixes : tuple, optional
        The prefixes to consider, defaults to ("http://", "https://")
    """

    def __init__(
        self,
        retries: int = 3,
        backoff_factor: float = 0.3,
        status_to_retry: Tuple[int, ...] = (500, 502, 504),
        prefixes: Tuple[str, ...] = ("http://", "https://"),
    ) -> None:
        self.session = Session()
        self.retries = retries

        r = Retry(
            total=retries,
            read=retries,
            connect=retries,
            backoff_factor=backoff_factor,
            status_forcelist=status_to_retry,
            method_whitelist=False,
        )
        adapter = HTTPAdapter(max_retries=r)
        for prefix in prefixes:
            self.session.mount(prefix, adapter)
        self.session.hooks = {"response": [lambda r, *args, **kwargs: r.raise_for_status()]}

    def get(self, url: str, payload: Optional[Mapping[str, Any]] = None,) -> Response:
        """Retrieve data from a url by GET and return the Response."""
        try:
            return self.session.get(url, params=payload)
        except (ConnectionError, HTTPError, RequestException, RetryError, Timeout):
            raise ConnectionError(f"Connection failed after {self.retries} retries.")
