import io
import shutil

import numpy as np
import pytest
import rasterio
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


@pytest.mark.flaky(max_runs=3)
def test_getmap_2file(geometry):
    dem = py3dep.get_map("DEM", geometry, 1e3, geo_crs=DEF_CRS, crs=ALT_CRS, output_dir="raster")
    with rasterio.open("raster/3DEPElevation:None_dd_0_0.gtiff") as f:
        mean = f.read().mean()

    shutil.rmtree("raster")
    assert abs(dem.mean().item() - mean) < 7e-2


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


def test_show_versions():
    f = io.StringIO()
    py3dep.show_versions(file=f)
    assert "INSTALLED VERSIONS" in f.getvalue()
