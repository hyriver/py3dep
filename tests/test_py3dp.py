import pytest
from shapely.geometry import Polygon

import py3dep


@pytest.fixture
def geom():
    return Polygon(
        [[-69.77, 45.07], [-69.31, 45.07], [-69.31, 45.45], [-69.77, 45.45], [-69.77, 45.07]]
    )


def test_dem(geom):
    dem = py3dep.get_map("DEM", geom, 10, geo_crs="epsg:4326", crs="epsg:3857", data_dir="data")
    dem = py3dep.get_map("DEM", geom.bounds, 1e3, geo_crs="epsg:4326", crs="epsg:3857")
    dem = py3dep.get_map("DEM", geom, 1e3, geo_crs="epsg:4326", crs="epsg:3857")
    gridxy = [dem.x.values, dem.y.values]
    e = py3dep.elevation_bygrid(gridxy=gridxy, resolution=1e3, crs=dem.crs.to_string().lower())
    assert abs((dem - e).sum().item() - (-14.198)) < 1e-3


def test_loc():
    elev = py3dep.elevation_byloc((-7766049.664788851, 5691929.739021257), "EPSG:3857")
    assert abs(elev - 356.59) < 1e-3


def test_deg2mpm(geom):
    slop = py3dep.get_map("Slope Degrees", geom, 1e3)
    slop = py3dep.utils.deg2mpm(slop)
    assert abs(slop.mean().item() - 0.05) < 1e-3
