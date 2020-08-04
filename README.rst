.. image:: https://raw.githubusercontent.com/cheginit/hydrodata/develop/docs/_static/py3dep_logo.png
    :target: https://github.com/cheginit/py3dep
    :align: center

|

=========== ===========================================================================
Package     Description
=========== ===========================================================================
Hydrodata_  Access NWIS, HCDN 2009, NLCD, and SSEBop databases
PyGeoOGC_   Query data from any ArcGIS RESTful-, WMS-, and WFS-based services
PyGeoUtils_ Convert responses from PyGeoOGC's supported web services to datasets
PyNHD_      Access NLDI and WaterData web services for navigating the NHDPlus database
Py3DEP_     Access topographic data through the 3D Elevation Program (3DEP) web service
PyDaymet_   Access the Daymet database for daily climate data
=========== ===========================================================================

.. _Hydrodata: https://github.com/cheginit/hydrodata
.. _PyGeoOGC: https://github.com/cheginit/pygeoogc
.. _PyGeoUtils: https://github.com/cheginit/pygeoutils
.. _PyNHD: https://github.com/cheginit/pynhd
.. _Py3DEP: https://github.com/cheginit/py3dep
.. _PyDaymet: https://github.com/cheginit/pydaymet

Py3DEP: Topographic data through 3DEP
-------------------------------------

.. image:: https://img.shields.io/pypi/v/py3dep.svg
    :target: https://pypi.python.org/pypi/py3dep
    :alt: PyPi

.. image:: https://img.shields.io/conda/vn/conda-forge/py3dep.svg
    :target: https://anaconda.org/conda-forge/py3dep
    :alt: Conda Version

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

🚨 **This package is under heavy development and breaking changes are likely to happen.** 🚨

Features
--------

Py3DEP provides access to the `3DEP <https://www.usgs.gov/core-science-systems/ngp/3dep>`__
database which is a part of the
`National Map services <https://viewer.nationalmap.gov/services/>`__.
The 3DEP service has multi-resolution sources and depending on the user provided resolution,
the data is resampled on the server-side based on all the available data sources. Py3DEP returns
the requests as `xarray <https://xarray.pydata.org/en/stable>`__ dataset. The 3DEP includes
the following layers:

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

Moreover, Py3DEP offers some additional utilities:

- ``elevation_bygrid``: For getting elevations of all the grid points in a 2D grid.
- ``elevation_byloc``: For getting elevation of a single point which is based on the National
  Map's `Elevation Point Query Service <https://nationalmap.gov/epqs/>`__.
- ``deg2mpm``: For converting slope dataset from degree to meter per meter.

You can try using Py3DEP without installing it on you system by clicking on the binder badge
below the Py3DEP banner. A Jupyter notebook instance with the Hydrodata software stack
pre-installed will be launched in your web browser and you can start coding!

Moreover, requests for additional functionalities can be submitted via
`issue tracker <https://github.com/cheginit/py3dep/issues>`__.


Installation
------------

You can install Py3DEP using ``pip`` after installing ``libgdal`` on your system
(for example, in Ubuntu run ``sudo apt install libgdal-dev``):

.. code-block:: console

    $ pip install py3dep

Alternatively, Py3DEP can be installed from the ``conda-forge`` repository
using `Conda <https://docs.conda.io/en/latest/>`__:

.. code-block:: console

    $ conda install -c conda-forge py3dep

Quick start
-----------

Py3DEP accepts `Shapely <https://shapely.readthedocs.io/en/latest/manual.html>`__'s
Polygon or a bounding box (a tuple of length four) as an input geometry.
We can use Hydrodata to get a watershed's geometry, then use it to get DEM and slope data
in meters/meters from Py3DEP using ``get_map`` function.

The ``get_map`` has a ``resolution`` argument that sets the target resolution
in meters. Note that the highest available resolution throughout the CONUS is about 10 m,
though higher resolutions are available in limited parts of the US. Note that the input
geometry can be in any valid spatial reference (``geo_crs`` argument). The ``crs`` argument,
however, is limited to ``CRS:84``, ``EPSG:4326``, and ``EPSG:3857`` since 3DEP only supports
these spatial references.

.. code-block:: python

    import py3dep
    from hydrodata import NLDI

    geom = NLDI().getfeature_byid("nwissite", "USGS-01031500", basin=True).geometry[0]
    dem = py3dep.get_map("DEM", geom, resolution=30, geo_crs="epsg:4326", crs="epsg:3857")
    slope = py3dep.get_map("Slope Degrees", geom, resolution=30)
    slope = py3dep.deg2mpm(slope)

.. image:: https://raw.githubusercontent.com/cheginit/hydrodata/develop/docs/_static/example_plots_py3dep.png
    :target: https://raw.githubusercontent.com/cheginit/hydrodata/develop/docs/_static/example_plots_py3dep.png
    :align: center

We can get the elevation for a single point within the US:

.. code-block:: python

    elev = py3dep.elevation_byloc((-7766049.665, 5691929.739), "epsg:3857")

Additionally, we can get the elevations of set of x- and y- coordinates of a grid. For example,
let's get the minimum temperature data within the watershed from Daymet using Hydrodata then
add the elevation as a new variable to the dataset:

.. code-block:: python

    import hydrodata.datasets as hds
    import xarray as xr
    import numpy as np

    clm = hds.daymet_bygeom(geom, dates=("2005-01-01", "2005-01-31"), variables="tmin")
    gridxy = (clm.x.values, clm.y.values)
    elev = py3dep.elevation_bygrid(gridxy, clm.crs, clm.res[0] * 1000)
    clm = xr.merge([clm, elev], combine_attrs="override")
    clm["elevation"] = clm.elevation.where(~np.isnan(clm.isel(time=0).tmin), drop=True)


Contributing
------------

Contributions are very welcomed. Please read
`CONTRIBUTING.rst <https://github.com/cheginit/pygeoogc/blob/master/CONTRIBUTING.rst>`__
file for instructions.
