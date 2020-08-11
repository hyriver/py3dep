import io

import numpy as np
import pytest
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
def test_dem(geometry):
    lyr = "DEM"
    py3dep.get_map(lyr, geometry.bounds, 1e3, geo_crs=DEF_CRS, crs=ALT_CRS)
    dem_10 = py3dep.get_map(lyr, geometry, 10, geo_crs=DEF_CRS, crs=ALT_CRS)
    dem_1e3 = py3dep.get_map(lyr, geometry, 1e3, geo_crs=DEF_CRS, crs=ALT_CRS)
    assert abs(abs(dem_10.mean().item() - dem_1e3.mean().item()) - 0.062) < 1e-3


def test_loc():
    elev = py3dep.elevation_byloc((-7766049.664788851, 5691929.739021257), ALT_CRS)
    assert abs(elev - 356.59) < 1e-3


@pytest.mark.flaky(max_runs=3)
def test_deg2mpm(geometry):
    slope = py3dep.get_map("Slope Degrees", geometry, 1e3)
    slope = py3dep.deg2mpm(slope)
    assert abs(slope.mean().item() - 0.05) < 1e-3


def test_grid(geometry):
    geo_crs = DEF_CRS
    crs = "+proj=lcc +lat_1=25 +lat_2=60 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0 +ellps=WGS84 +units=km +no_defs"
    geom = MatchCRS.geometry(geometry, geo_crs, crs)
    xmin, ymin, xmax, ymax = geom.bounds
    res = 1
    gx = np.arange(xmin, xmax, res)
    gy = np.arange(ymin, ymax, res)
    elev = py3dep.elevation_bygrid((gx, gy), crs, res * 1e3)
    assert abs(elev.mean().item() - 295.686) < 1e-3


def test_show_versions():
    f = io.StringIO()
    py3dep.show_versions(file=f)
    assert "INSTALLED VERSIONS" in f.getvalue()
