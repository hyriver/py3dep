import io
import shutil
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import pyproj
import pytest
import rioxarray as rxr
from pygeoogc import utils
from shapely import ops
from shapely.geometry import MultiLineString, Polygon

import py3dep
from py3dep.cli import cli

DEF_CRS = "epsg:4326"
ALT_CRS = "epsg:3857"
GEOM = Polygon(
    [[-69.77, 45.07], [-69.31, 45.07], [-69.31, 45.45], [-69.77, 45.45], [-69.77, 45.07]]
)
LINE = MultiLineString(
    [
        [[-69.77, 45.07], [-69.31, 45.07]],
        [[-69.31, 45.07], [-69.31, 45.45]],
        [[-69.31, 45.45], [-69.77, 45.45]],
    ]
)
LYR = "Slope Degrees"
SMALL = 1e-3


def test_profile():
    ep = py3dep.elevation_profile(LINE, 10)
    epm = py3dep.elevation_profile(ops.linemerge(LINE), 10)
    assert (
        abs(ep.mean().compute().item() - 389.7205) < SMALL
        and abs(epm.mean().compute().item() - 389.7205) < SMALL
    )


def test_getmap():
    layers = ["DEM", LYR]
    ds = py3dep.get_map(layers, GEOM.bounds, 1e3, geo_crs=DEF_CRS, crs=ALT_CRS)
    dem_10 = py3dep.get_map(layers[0], GEOM, 10, geo_crs=DEF_CRS, crs=ALT_CRS)
    dem_1e3 = py3dep.get_map(layers[0], GEOM, 1e3, geo_crs=DEF_CRS, crs=ALT_CRS)
    fpath = Path("dem_10.tif")
    dem_10.rio.to_raster(fpath)
    dem_10 = rxr.open_rasterio(fpath)
    assert (
        sorted(ds.keys()) == ["elevation", "slope_degrees"]
        and abs(dem_10.mean().compute().item() - dem_1e3.mean().compute().item()) < 7e-2
    )
    dem_10.close()
    fpath.unlink()


def test_fill_depressions():
    ds = py3dep.get_map("DEM", GEOM.bounds, 1e3)
    ds = py3dep.fill_depressions(ds)
    assert abs(ds.mean().compute().item() - 296.206) < SMALL


@pytest.mark.parametrize(
    "source,expected",
    [("airmap", 363), ("tnm", 356.59), ("tep", 356.088)],
)
def test_bycoords(source, expected):
    coords = [(-7766049.664788851, 5691929.739021257)]
    dem = py3dep.elevation_bycoords(coords * 101, pyproj.CRS(ALT_CRS), source=source)
    assert abs(sum(set(dem)) - expected) < SMALL


def test_deg2mpm():
    slope = py3dep.get_map(LYR, GEOM, 1e3)
    slope = py3dep.deg2mpm(slope)
    assert abs(slope.mean().compute().item() - 0.050) < SMALL


def test_grid():
    crs = "+proj=lcc +lat_1=25 +lat_2=60 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"
    geom = utils.match_crs(GEOM, DEF_CRS, crs)
    xmin, ymin, xmax, ymax = geom.bounds
    res = 1e3
    gx = np.arange(xmin, xmax, res)
    gy = np.arange(ymin, ymax, res)
    elev = py3dep.elevation_bygrid(tuple(gx), tuple(gy), crs, res)
    elev_fill = py3dep.elevation_bygrid(tuple(gx), tuple(gy), crs, res, depression_filling=True)
    assert ((elev_fill - elev).sum().compute().item() - 2710.839) < SMALL


def test_check_3dep_availability():
    avail = py3dep.check_3dep_availability(GEOM.bounds)
    assert avail["1m"] and avail["10m"] and avail["30m"]


def test_query_3dep_source():
    src = py3dep.query_3dep_sources(GEOM.bounds)
    res_all = src.groupby("dem_res")["OBJECTID"].count().to_dict()
    src = py3dep.query_3dep_sources(GEOM.bounds, res="1m")
    res_1m = src.groupby("dem_res")["OBJECTID"].count().to_dict()
    assert res_all == {"10m": 8, "1m": 4, "30m": 8} and res_1m == {"1m": 4}


class TestCLI:
    """Test the command-line interface."""

    def test_geometry(self, runner):
        gdf = gpd.GeoDataFrame(
            {"id": "geo_test", "res": 1e3}, geometry=[GEOM], index=[0], crs=DEF_CRS
        )
        geo_gpkg = Path("nat_geo.gpkg")
        gdf.to_file(geo_gpkg)
        ret = runner.invoke(
            cli, ["geometry", str(geo_gpkg), "-l", "DEM", "-l", LYR, "-s", "geo_map"]
        )
        if geo_gpkg.is_dir():
            shutil.rmtree(geo_gpkg)
        else:
            geo_gpkg.unlink()
        assert ret.exit_code == 0
        assert "Found 1 geometry" in ret.output
        shutil.rmtree("geo_map")

    def test_coords(self, runner):
        df = pd.DataFrame(
            [[-69.77, 45.07], [-69.31, 45.07], [-69.31, 45.45]], columns=["lon", "lat"]
        )
        coord_csv = "coords.csv"
        df.to_csv(coord_csv)
        ret = runner.invoke(cli, ["coords", coord_csv, "-s", "geo_coords", "-q", "airmap"])
        Path(coord_csv).unlink()
        assert ret.exit_code == 0
        assert "Found coordinates of 3 points" in ret.output
        shutil.rmtree("geo_coords")


def test_show_versions():
    f = io.StringIO()
    py3dep.show_versions(file=f)
    assert "INSTALLED VERSIONS" in f.getvalue()
