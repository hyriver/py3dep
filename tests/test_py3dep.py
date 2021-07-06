import io
import shutil
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from pygeoogc import utils
from shapely.geometry import Polygon

import py3dep

DEF_CRS = "epsg:4326"
ALT_CRS = "epsg:3857"
GEOM = Polygon(
    [[-69.77, 45.07], [-69.31, 45.07], [-69.31, 45.45], [-69.77, 45.45], [-69.77, 45.07]]
)
LYR = "Slope Degrees"


def test_getmap():
    layers = ["DEM", LYR]
    ds = py3dep.get_map(layers, GEOM.bounds, 1e3, geo_crs=DEF_CRS, crs=ALT_CRS)
    dem_10 = py3dep.get_map(layers[0], GEOM, 10, geo_crs=DEF_CRS, crs=ALT_CRS)
    dem_1e3 = py3dep.get_map(layers[0], GEOM, 1e3, geo_crs=DEF_CRS, crs=ALT_CRS)
    assert (
        sorted(ds.keys()) == ["elevation", "slope_degrees"]
        and abs(dem_10.mean().compute().item() - dem_1e3.mean().compute().item()) < 7e-2
    )


def test_coords():
    airmap = py3dep.elevation_bycoords([(-7766049.664788851, 5691929.739021257)] * 200, ALT_CRS)
    tnm = py3dep.elevation_bycoords(
        [(-7766049.664788851, 5691929.739021257)] * 200, ALT_CRS, source="tnm"
    )
    assert set(tnm) == {1169.9} and set(airmap) == {363}


def test_deg2mpm():
    slope = py3dep.get_map(LYR, GEOM, 1e3)
    slope = py3dep.deg2mpm(slope)
    assert abs(slope.mean().compute().item() - 0.05) < 1e-3


def test_grid():
    geo_crs = DEF_CRS
    crs = "+proj=lcc +lat_1=25 +lat_2=60 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"
    geom = utils.match_crs(GEOM, geo_crs, crs)
    xmin, ymin, xmax, ymax = geom.bounds
    res = 1e3
    gx = np.arange(xmin, xmax, res)
    gy = np.arange(ymin, ymax, res)
    elev = py3dep.elevation_bygrid(gx, gy, crs, res)
    elev_fill = py3dep.elevation_bygrid(gx, gy, crs, res, depression_filling=True)
    assert ((elev_fill - elev).sum().item() - 1404.618) < 1e-3


def test_cli_map(script_runner):
    gdf = gpd.GeoDataFrame({"id": "geo_test", "res": 1e3}, geometry=[GEOM], index=[0])
    geo_gpkg = "nat_geo.gpkg"
    gdf.to_file(geo_gpkg)
    ret = script_runner.run("py3dep", geo_gpkg, "geometry", DEF_CRS, "-l", LYR, "-s", "geo_map")
    shutil.rmtree(geo_gpkg)
    shutil.rmtree("geo_map")
    assert ret.success
    assert "Retrieved topography data for 1 item(s)." in ret.stdout
    assert ret.stderr == ""


def test_cli_coords(script_runner):
    df = pd.DataFrame([(-7766049.664788851, 5691929.739021257)] * 3, columns=["x", "y"])
    coord_csv = "coords.csv"
    df.to_csv(coord_csv)
    ret = script_runner.run("py3dep", coord_csv, "coords", ALT_CRS, "-s", "geo_coords")
    Path(coord_csv).unlink()
    shutil.rmtree("geo_coords")
    assert ret.success
    assert "Retrieved elevation data for 3 item(s)." in ret.stdout
    assert ret.stderr == ""


def test_show_versions():
    f = io.StringIO()
    py3dep.show_versions(file=f)
    assert "INSTALLED VERSIONS" in f.getvalue()
