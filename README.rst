.. image:: https://raw.githubusercontent.com/cheginit/hydrodata/develop/docs/_static/py3dep_logo.png
    :target: https://github.com/cheginit/py3dep
    :align: center

|

.. image:: https://img.shields.io/pypi/v/py3dep.svg
    :target: https://pypi.python.org/pypi/py3dep
    :alt: PyPi

.. image:: https://codecov.io/gh/cheginit/py3dep/branch/master/graph/badge.svg
    :target: https://codecov.io/gh/cheginit/py3dep
    :alt: CodeCov

.. image:: https://github.com/cheginit/py3dep/workflows/build/badge.svg
    :target: https://github.com/cheginit/py3dep/workflows/build
    :alt: Github Actions

.. image:: https://mybinder.org/badge_logo.svg
    :target: https://mybinder.org/v2/gh/cheginit/hydrodata/develop
    :alt: Binder

|

.. image:: https://www.codefactor.io/repository/github/cheginit/py3dep/badge
   :target: https://www.codefactor.io/repository/github/cheginit/py3dep
   :alt: CodeFactor

.. image:: https://img.shields.io/badge/code%20style-black-000000.svg
    :target: https://github.com/psf/black
    :alt: black

.. image:: https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white
    :target: https://github.com/pre-commit/pre-commit
    :alt: pre-commit

|

Features
--------

Py3DEP is a part of `Hydrodata <https://github.com/cheginit/hydrodata>`__ software stack
and provides access to `3DEP <https://www.usgs.gov/core-science-systems/ngp/3dep>`__ which is
a part the National Map services. The 3DEP service has multi-resolution sources
and depending on the user provided resolution, the data is resampled on
server-side based on all the available data sources. Py3DEP returns All the requestes
as ``xarray.DataArray`` or ``xarray.Dataset`` that allows
the user to process the data using all the ``xarray``'s functionalities.
The following layers are available:

- DEM
- Hillshade Gray
- Aspect Degrees
- Aspect Map
- GreyHillshade Elevation Fill
- Hillshade Multidirectional
- Slope Map
- Slope Degrees
- Hillshade Elevation Tinted
- Height Ellipsoidal
- Contour 25
- Contour Smoothed 25

Moreover, Py3DEP offers thress additonal function:

- ``elevation_bygrid``: For getting elevations of all the grid points in a 2D array of
  x- and y-coordinates.
- ``elevation_byloc``: For getting elevation of a single point based on the National
  Map's `Elevation Point Query Service <https://nationalmap.gov/epqs/>`__ service.
- ``deg2mpm``: For converting slop from degree to meter per meter.

Moreover, requests for additional databases or functionalities can be submitted via
`issue tracker <https://github.com/cheginit/py3dep/issues>`__.


Installation
------------

You can install py3dep using ``pip``:

.. code-block:: console

    $ pip install py3dep

Quickstart
----------

Py3DEP accepts `Shapely <https://shapely.readthedocs.io/en/latest/manual.html>`__'s
Polygon or a bounding box (a tuple of length four) as input an input geometry.
We can use `Hydrodata <https://hydrodata.readthedocs.io/en/latest/>`__ to get
a watershed's geometry, then use this geometry to get DEM and slope data
in meters/meters from Py3DEP. The ``get_map`` function has a ``resolution`` argument
that sets the target resolution in meters. Note that the highest available resolution
throughout the CONUS is about 10 m, though higher resolutions are available for limited
parts of the US. The geometry can be in any valid spatial reference (``geo_crs``).
The ``crs`` argument, however, is limited to ``CRS:84``, ``EPSG:4326``, and ``EPSG:3857``
since 3DEP only supports these spatial references.

.. code-block:: python

    import py3dep
    from hydrodata import NLDI

    geom = NLDI().getfeature_byid("nwissite", "USGS-01031500", basin=True).geometry[0]
    dem = py3dep.get_map("DEM", geom, resolution=30, geo_crs="epsg:4326", crs="epsg:3857")
    slope = py3dep.get_map("Slope Degrees", geom, resolution=30)
    slope = py3dep.utils.deg2mpm(slope)

.. image:: https://raw.githubusercontent.com/cheginit/py3dep/master/docs/_static/example_plots.png
    :target: https://raw.githubusercontent.com/cheginit/py3dep/master/docs/_static/example_plots.png
    :align: center

We can get the elevation of a single point within the US:

.. code-block:: python

    elev = py3dep.elevation_byloc((-7766049.665, 5691929.739), "epsg:3857")

Additionally, we can get the elevations of set of x- and y- coordinates of a grid. For example,
let's get the minimum temperature data within the watershed using Hydrodata then add the elevation
as a new variable:

.. code-block:: python

    import hydrodata.datasets as hds
    import xarray as xr
    import numpy as np

    clm = hds.daymet_bygeom(geom, dates=("2005-01-01", "2005-01-31"), variables="tmin")
    gridxy = (clm.x.values, clm.y.values)
    res = clm.res[0] * 1000
    elev = py3dep.elevation_bygrid(gridxy, clm.crs, res)
    clm = xr.merge([clm, elev], combine_attrs="override")
    clm["elevation"] = clm.elevation.where(~np.isnan(clm.isel(time=0).tmin), drop=True)


Contributing
------------

If you are interested in contributing to the project please get in touch.
You can find information about contributing to py3dep at our
`Contributing page <https://py3dep.readthedocs.io/en/latest/contributing.html>`__.

