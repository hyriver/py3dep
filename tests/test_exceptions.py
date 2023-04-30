import sys

import pytest
from shapely.geometry import Polygon

import py3dep
from py3dep import InputTypeError, InputValueError

DEF_CRS = 4326
GEOM = Polygon(
    [[-69.77, 45.07], [-69.31, 45.07], [-69.31, 45.45], [-69.77, 45.45], [-69.77, 45.07]]
)


has_typeguard = True if sys.modules.get("typeguard") else False


@pytest.mark.skipif(has_typeguard, reason="Broken if Typeguard is enabled")
def test_wrong_bbox():
    with pytest.raises(InputTypeError) as ex:
        _ = py3dep.check_3dep_availability((1, 2, 3))
        assert "tuple of length 4" in str(ex.value)


def test_wrong_layer():
    with pytest.raises(InputValueError) as ex:
        _ = py3dep.get_map("None", GEOM.bounds, 1e3)
        assert "DEM" in str(ex.value)


def test_wrong_crs():
    with pytest.raises(InputValueError) as ex:
        _ = py3dep.get_map("None", GEOM.bounds, 1e3, DEF_CRS, "ESRI:102003")
        assert "crs" in str(ex.value)
