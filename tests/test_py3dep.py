import io
import shutil
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
import rasterio
import xarray as xr
from pygeoogc import MatchCRS
from shapely.geometry import Polygon

import py3dep

DEF_CRS = "epsg:4326"
ALT_CRS = "epsg:3857"


@pytest.fixture
def geometry():
    return Polygon(
        [[-69.77, 45.07], [-69.31, 45.07], [-69.31, 45.45], [-69.77, 45.45], [-69.77, 45.07]]
    )


@pytest.mark.flaky(max_runs=3)
def test_getmap(geometry):
    lyr = ["DEM", "Slope Degrees"]
    ds = py3dep.get_map(lyr, geometry.bounds, 1e3, geo_crs=DEF_CRS, crs=ALT_CRS)
    dem_10 = py3dep.get_map(lyr[0], geometry, 10, geo_crs=DEF_CRS, crs=ALT_CRS)
    dem_1e3 = py3dep.get_map(lyr[0], geometry, 1e3, geo_crs=DEF_CRS, crs=ALT_CRS)
    assert (
        sorted(ds.keys()) == ["elevation", "slope_degrees"]
        and abs(dem_10.mean().item() - dem_1e3.mean().item()) < 7e-2
    )


def test_coords():
    elev = py3dep.elevation_bycoords([(-7766049.664788851, 5691929.739021257)] * 3, ALT_CRS)
    assert elev == [363] * 3


@pytest.mark.flaky(max_runs=3)
def test_deg2mpm(geometry):
    slope = py3dep.get_map("Slope Degrees", geometry, 1e3)
    slope = py3dep.deg2mpm(slope)
    assert abs(slope.mean().item() - 0.05) < 1e-3


def test_grid(geometry):
    geo_crs = DEF_CRS
    crs = "+proj=lcc +lat_1=25 +lat_2=60 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"
    geom = MatchCRS.geometry(geometry, geo_crs, crs)
    xmin, ymin, xmax, ymax = geom.bounds
    res = 1e3
    gx = np.arange(xmin, xmax, res)
    gy = np.arange(ymin, ymax, res)
    elev = py3dep.elevation_bygrid(gx, gy, crs, res)
    assert abs(elev.mean().item() - 295.763) < 1e-3


def test_cli_map(script_runner, geometry):
    gdf = gpd.GeoDataFrame({"id": "geo_test", "res": 1e3}, geometry=[geometry], index=[0])
    gdf.to_file("nat_geo.gpkg")
    ret = script_runner.run(
        "py3dep", "nat_geo.gpkg", "geometry", DEF_CRS, "-l", "Slope Degrees", "-s", "geo_map"
    )
    shutil.rmtree("nat_geo.gpkg")
    shutil.rmtree("geo_map")
    assert ret.success
    assert "Retrieved topography data for 1 item(s)." in ret.stdout
    assert ret.stderr == ""


def test_cli_coords(script_runner):
    df = pd.DataFrame([(-7766049.664788851, 5691929.739021257)] * 3, columns=["x", "y"])
    df.to_csv("coords.csv")
    ret = script_runner.run("py3dep", "coords.csv", "coords", ALT_CRS, "-s", "geo_coords")
    Path("coords.csv").unlink()
    shutil.rmtree("geo_coords")
    assert ret.success
    assert "Retrieved elevation data for 3 item(s)." in ret.stdout
    assert ret.stderr == ""


def test_show_versions():
    f = io.StringIO()
    py3dep.show_versions(file=f)
    assert "INSTALLED VERSIONS" in f.getvalue()
