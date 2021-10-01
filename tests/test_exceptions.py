import numpy as np
import pytest
import xarray as xr

import py3dep
from py3dep import MissingAttribute, MissingColumns, MissingCRS, MissingDependency


def test_missing_nodata():
    with pytest.raises(MissingAttribute) as ex:
        dem = xr.DataArray(np.random.randn(2, 3), dims=("x", "y"), coords={"x": [10, 20]})
        _ = py3dep.fill_depressions(dem)
        assert "nodatavals" in str(ex.value)
