import pytest
from shapely.geometry import Polygon

import py3dep
from py3dep import InvalidInputType, InvalidInputValue

DEF_CRS = "epsg:4326"
GEOM = Polygon(
    [[-69.77, 45.07], [-69.31, 45.07], [-69.31, 45.45], [-69.77, 45.45], [-69.77, 45.07]]
)

try:
    import typeguard  # noqa: F401
except ImportError:
    has_typeguard = False
else:
    has_typeguard = True


@pytest.mark.skipif(has_typeguard, reason="Broken if Typeguard is enabled")
def test_wrong_bbox():
    with pytest.raises(InvalidInputType) as ex:
        _ = py3dep.check_3dep_availability((1, 2, 3))
        assert "tuple of length 4" in str(ex.value)


def test_wrong_layer():
    with pytest.raises(InvalidInputValue) as ex:
        _ = py3dep.get_map("None", GEOM.bounds, 1e3)
        assert "DEM" in str(ex.value)


def test_wrong_crs():
    with pytest.raises(InvalidInputValue) as ex:
        _ = py3dep.get_map("None", GEOM.bounds, 1e3, DEF_CRS, "ESRI:102003")
        assert "crs" in str(ex.value)
