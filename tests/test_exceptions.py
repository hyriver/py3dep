import pytest
from shapely import Polygon

import py3dep
from py3dep import InputTypeError, InputValueError

DEF_CRS = 4326
GEOM = Polygon(
    [[-69.77, 45.07], [-69.31, 45.07], [-69.31, 45.45], [-69.77, 45.45], [-69.77, 45.07]]
)


def test_wrong_bbox():
    bbox = (1, 2, 3)
    with pytest.raises(InputTypeError, match="tuple of length 4"):
        _ = py3dep.check_3dep_availability(bbox)


def test_wrong_layer():
    geom = GEOM.bounds
    with pytest.raises(InputValueError, match="DEM"):
        _ = py3dep.get_map("None", geom, 1e3)


def test_wrong_crs():
    geom = GEOM.bounds
    with pytest.raises(InputValueError, match="crs"):
        _ = py3dep.get_map("DEM", geom, 1e3, DEF_CRS, "ESRI:102003")
