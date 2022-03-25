"""Utilities for Py3DEP."""
from typing import List, Tuple, TypeVar, Union

import geopandas as gpd
import numpy as np
import pygeoutils as geoutils
import rioxarray  # noqa: F401
import xarray as xr
from pygeoogc import ArcGISRESTful
from shapely.geometry import Polygon

try:
    import richdem as rd
except ImportError:
    rd = None

from .exceptions import MissingDependency

__all__ = ["deg2mpm", "fill_depressions"]
X = TypeVar("X", xr.DataArray, xr.Dataset)
DEF_CRS = "EPSG:4326"


def fill_depressions(dem_da: xr.DataArray) -> xr.DataArray:
    """Fill depressions and adjust flat areas in a DEM using `RichDEM <https://richdem.readthedocs.io>`__.

    Parameters
    ----------
    dem : xarray.DataArray or numpy.ndarray
        Digital Elevation Model.

    Returns
    -------
    xarray.DataArray
        Conditioned DEM after applying
        `depression filling <https://richdem.readthedocs.io/en/latest/depression_filling.html>`__
        and
        `flat area resolution <https://richdem.readthedocs.io/en/latest/flat_resolution.html>`__
        operations.
    """
    dem = dem_da.copy()
    if rd is None:
        raise MissingDependency

    nodata = -9999
    with xr.set_options(keep_attrs=True):  # type: ignore
        dem = dem.fillna(nodata)
        rda = rd.rdarray(dem, no_data=nodata)
        rda.projection = dem.rio.crs
        rda.geotransform = geoutils.transform2tuple(dem.rio.transform())
        rda = rd.FillDepressions(rda, epsilon=False)
        dem.data = rd.ResolveFlats(rda)
        return dem.where(dem > nodata, dem.attrs["nodatavals"][0])  # type: ignore


def deg2mpm(slope: xr.DataArray) -> xr.DataArray:
    """Convert slope from degrees to meter/meter.

    Parameters
    ----------
    slope : xarray.DataArray
        Slope in degrees.

    Returns
    -------
    xarray.DataArray
        Slope in meter/meter. The name is set to ``slope`` and the ``units`` attribute
        is set to ``m/m``.
    """
    with xr.set_options(keep_attrs=True):  # type: ignore
        if hasattr(slope, "_FillValue"):
            nodata = slope.attrs["_FillValue"]
        else:
            nodata = slope.attrs["nodatavals"][0]
        slope = slope.where(slope != nodata, drop=False)
        slope = np.tan(np.deg2rad(slope))
        slope.attrs["nodatavals"] = (np.nan,)
        if hasattr(slope, "_FillValue"):
            slope.attrs["_FillValue"] = np.nan
        slope.name = "slope"
        slope.attrs["units"] = "m/m"
    return slope


def rename_layers(ds: X, valid_layers: List[str]) -> X:
    """Rename layers in a dataset."""
    rename = {lyr: lyr.split(":")[-1].replace(" ", "_").lower() for lyr in valid_layers}
    rename.update({"3DEPElevation:None": "elevation"})

    if isinstance(ds, xr.DataArray):
        ds.name = rename[str(ds.name)]
    else:
        ds = ds.rename({n: rename[str(n)] for n in ds})

    return ds


class RESTful:
    """Base class for getting geospatial data from a ArcGISRESTful service.

    Parameters
    ----------
    base_url : str, optional
        The ArcGIS RESTful service url. The URL must either include a layer number
        after the last ``/`` in the url or the target layer must be passed as an argument.
    layer : int, optional
        A valid service layer.
    outfields : str or list, optional
        Target field name(s), default to "*" i.e., all the fields.
    crs : str, optional
        Target spatial reference, default to EPSG:4326
    outformat : str, optional
        One of the output formats offered by the selected layer. If not correct
        a list of available formats is shown, defaults to ``json``.
    """

    def __init__(
        self,
        base_url: str,
        layer: int,
        outfields: Union[str, List[str]] = "*",
        crs: str = DEF_CRS,
        outformat: str = "json",
    ) -> None:
        self.client = ArcGISRESTful(
            base_url,
            layer,
            outformat=outformat,
            outfields=outfields,
            crs=crs,
        )

    def bygeom(
        self,
        geom: Union[Polygon, List[Tuple[float, float]], Tuple[float, float, float, float]],
        geo_crs: str = DEF_CRS,
    ) -> gpd.GeoDataFrame:
        """Get feature within a geometry that can be combined with a SQL where clause.

        Parameters
        ----------
        geom : Polygon or tuple
            A geometry (Polygon) or bounding box (tuple of length 4).
        geo_crs : str
            The spatial reference of the input geometry.

        Returns
        -------
        geopandas.GeoDataFrame
            The requested features as a GeoDataFrame.
        """
        oids = self.client.oids_bygeom(geom, geo_crs=geo_crs)
        return geoutils.json2geodf(self.client.get_features(oids), self.client.client.crs)

    def bysql(
        self,
        sql_clause: str,
    ) -> gpd.GeoDataFrame:
        """Get feature IDs using a valid SQL 92 WHERE clause.

        Notes
        -----
        Not all web services support this type of query. For more details look
        `here <https://developers.arcgis.com/rest/services-reference/query-feature-service-.htm#ESRI_SECTION2_07DD2C5127674F6A814CE6C07D39AD46>`__

        Parameters
        ----------
        sql_clause : str
            A valid SQL 92 WHERE clause.

        Returns
        -------
        geopandas.GeoDataFrame
            The requested features as a GeoDataFrame.
        """
        oids = self.client.oids_bysql(sql_clause)
        return geoutils.json2geodf(self.client.get_features(oids), self.client.client.crs)
