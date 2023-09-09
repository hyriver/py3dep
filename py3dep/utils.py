"""Utilities for Py3DEP."""
from __future__ import annotations

from typing import TYPE_CHECKING, Sequence, Union, overload

import geopandas as gpd
import numpy as np
import pygeoutils as geoutils
import pyproj
import xarray as xr
from pygeoogc import ArcGISRESTful

try:
    import pyflwdir
except ImportError:
    pyflwdir = None

if TYPE_CHECKING:
    from shapely.geometry import Polygon

    CRSTYPE = Union[int, str, pyproj.CRS]

from py3dep.exceptions import DependencyError

__all__ = ["deg2mpm", "fill_depressions"]


def fill_depressions(dem: xr.DataArray) -> xr.DataArray:
    """Fill depressions and adjust flat areas in a DEM using `RichDEM <https://richdem.readthedocs.io>`__.

    Parameters
    ----------
    dem_da : xarray.DataArray or numpy.ndarray
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
    dem = dem.astype("f8")
    if pyflwdir is None:
        raise DependencyError

    attrs = dem.attrs
    filled, _ = pyflwdir.dem.fill_depressions(dem.values, outlets="min", nodata=np.nan)
    dem = dem.copy(data=filled)
    dem.attrs = attrs
    dem = dem.rio.write_nodata(np.nan)
    return dem


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
    with xr.set_options(keep_attrs=True):
        if hasattr(slope, "_FillValue"):
            nodata = slope.attrs["_FillValue"]
        elif hasattr(slope, "nodatavals"):
            _nodata = slope.attrs["nodatavals"]
            nodata = _nodata[0] if isinstance(_nodata, Sequence) else _nodata
        else:
            nodata = np.nan
        slope = slope.where(slope != nodata, drop=False)

        def to_mpm(da: xr.DataArray) -> xr.DataArray:
            """Convert slope from degrees to meter/meter."""
            return xr.apply_ufunc(
                lambda x: np.tan(np.deg2rad(x)),
                da,
                vectorize=True,
            )

        slope = to_mpm(slope).compute()
        slope.attrs["nodatavals"] = (np.nan,)
        if hasattr(slope, "_FillValue"):
            slope.attrs["_FillValue"] = np.nan
        slope.name = "slope"
        slope.attrs["units"] = "m/m"
    return slope


@overload
def rename_layers(ds: xr.DataArray, valid_layers: list[str]) -> xr.DataArray:
    ...


@overload
def rename_layers(ds: xr.Dataset, valid_layers: list[str]) -> xr.Dataset:
    ...


def rename_layers(
    ds: xr.DataArray | xr.Dataset, valid_layers: list[str]
) -> xr.DataArray | xr.Dataset:
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
    crs : str, int, or pyproj.CRS, optional
        Target spatial reference, default to EPSG:4326
    outformat : str, optional
        One of the output formats offered by the selected layer. If not correct
        a list of available formats is shown, defaults to ``json``.
    """

    def __init__(
        self,
        base_url: str,
        layer: int,
        outfields: str | list[str] = "*",
        crs: CRSTYPE = 4326,
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
        geom: Polygon | list[tuple[float, float]] | tuple[float, float, float, float],
        geo_crs: CRSTYPE = 4326,
    ) -> gpd.GeoDataFrame:
        """Get feature within a geometry that can be combined with a SQL where clause.

        Parameters
        ----------
        geom : Polygon or tuple
            A geometry (Polygon) or bounding box (tuple of length 4).
        geo_crs : str, int, or pyproj.CRS
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
