.. image:: https://raw.githubusercontent.com/cheginit/HyRiver-examples/main/notebooks/_static/py3dep_logo.png
    :target: https://github.com/cheginit/HyRiver

|

.. |pygeohydro| image:: https://github.com/cheginit/pygeohydro/actions/workflows/test.yml/badge.svg
    :target: https://github.com/cheginit/pygeohydro/actions?query=workflow%3Apytest
    :alt: Github Actions

.. |pygeoogc| image:: https://github.com/cheginit/pygeoogc/actions/workflows/test.yml/badge.svg
    :target: https://github.com/cheginit/pygeoogc/actions?query=workflow%3Apytest
    :alt: Github Actions

.. |pygeoutils| image:: https://github.com/cheginit/pygeoutils/actions/workflows/test.yml/badge.svg
    :target: https://github.com/cheginit/pygeoutils/actions?query=workflow%3Apytest
    :alt: Github Actions

.. |pynhd| image:: https://github.com/cheginit/pynhd/actions/workflows/test.yml/badge.svg
    :target: https://github.com/cheginit/pynhd/actions?query=workflow%3Apytest
    :alt: Github Actions

.. |py3dep| image:: https://github.com/cheginit/py3dep/actions/workflows/test.yml/badge.svg
    :target: https://github.com/cheginit/py3dep/actions?query=workflow%3Apytest
    :alt: Github Actions

.. |pydaymet| image:: https://github.com/cheginit/pydaymet/actions/workflows/test.yml/badge.svg
    :target: https://github.com/cheginit/pydaymet/actions?query=workflow%3Apytest
    :alt: Github Actions

=========== ==================================================================== ============
Package     Description                                                          Status
=========== ==================================================================== ============
PyGeoHydro_ Access NWIS, NID, HCDN 2009, NLCD, and SSEBop databases              |pygeohydro|
PyGeoOGC_   Send queries to any ArcGIS RESTful-, WMS-, and WFS-based services    |pygeoogc|
PyGeoUtils_ Convert responses from PyGeoOGC's supported web services to datasets |pygeoutils|
PyNHD_      Navigate and subset NHDPlus (MR and HR) using web services           |pynhd|
Py3DEP_     Access topographic data through National Map's 3DEP web service      |py3dep|
PyDaymet_   Access Daymet for daily climate data both single pixel and gridded   |pydaymet|
=========== ==================================================================== ============

.. _PyGeoHydro: https://github.com/cheginit/pygeohydro
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

.. image:: https://mybinder.org/badge_logo.svg
    :target: https://mybinder.org/v2/gh/cheginit/pygeohydro/master?filepath=docs%2Fexamples
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

Py3DEP is a part of `HyRiver <https://github.com/cheginit/HyRiver>`__ software stack that
is designed to aid in watershed analysis through web services. This package provides
access to the `3DEP <https://www.usgs.gov/core-science-systems/ngp/3dep>`__
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
- ``elevation_bycoords``: For getting elevation of a list of ``x`` and ``y`` coordinates.
- ``deg2mpm``: For converting slope dataset from degree to meter per meter.

You can try using Py3DEP without installing it on you system by clicking on the binder badge
below the Py3DEP banner. A Jupyter notebook instance with the stack
pre-installed will be launched in your web browser and you can start coding!

Please note that since this project is in early development stages, while the provided
functionalities should be stable, changes in APIs are possible in new releases. But we
appreciate it if you give this project a try and provide feedback. Contributions are most welcome.

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
We can use PyNHD to get a watershed's geometry, then use it to get the DEM and slope
in meters/meters from Py3DEP using ``get_map`` function.

The ``get_map`` has a ``resolution`` argument that sets the target resolution
in meters. Note that the highest available resolution throughout the CONUS is about 10 m,
though higher resolutions are available in limited parts of the US. Note that the input
geometry can be in any valid spatial reference (``geo_crs`` argument). The ``crs`` argument,
however, is limited to ``CRS:84``, ``EPSG:4326``, and ``EPSG:3857`` since 3DEP only supports
these spatial references.

.. code-block:: python

    import py3dep
    from pynhd import NLDI

    geom = NLDI().get_basins("01031500").geometry[0]
    dem = py3dep.get_map("DEM", geom, resolution=30, geo_crs="epsg:4326", crs="epsg:3857")
    slope = py3dep.get_map("Slope Degrees", geom, resolution=30)
    slope = py3dep.deg2mpm(slope)

.. image:: https://raw.githubusercontent.com/cheginit/HyRiver-examples/main/notebooks/_static/dem_slope.png
    :target: https://github.com/cheginit/HyRiver-examples/blob/main/notebooks/3dep.ipynb
    :align: center

The ``get_map`` function also has another argument for saving the dataset into a raster file. We
should provide the path to a folder:

.. code-block:: python

    dem = py3dep.get_map("DEM", geom, 1e3, output_dir="raster")

Moreover, we can get the elevations of set of x- and y- coordinates on a grid. For example,
let's get the minimum temperature data within this watershed from Daymet using PyDaymet then
add the elevation as a new variable to the dataset:

.. code-block:: python

    import pydaymet as daymet
    import xarray as xr
    import numpy as np

    clm = daymet.get_bygeom(geometry, ("2005-01-01", "2005-01-31"), variables="tmin")
    elev = py3dep.elevation_bygrid(clm.x.values, clm.y.values, clm.crs, clm.res[0] * 1000)
    attrs = clm.attrs
    clm = xr.merge([clm, elev])
    clm["elevation"] = clm.elevation.where(~np.isnan(clm.isel(time=0).tmin), drop=True)
    clm.attrs.update(attrs)

Now, let's get street network data using `osmnx <https://github.com/gboeing/osmnx>`_ package
and add elevation data for its nodes using ``elevation_bycoords`` function.

.. code-block:: python

    import osmnx as ox

    G = ox.graph_from_place("Piedmont, California, USA", network_type="drive")
    x, y = nx.get_node_attributes(G, "x").values(), nx.get_node_attributes(G, "y").values()
    elevation = py3dep.elevation_bycoords(zip(x, y), crs="epsg:4326")
    nx.set_node_attributes(G, dict(zip(G.nodes(), elevation)), "elevation")

.. image:: https://raw.githubusercontent.com/cheginit/HyRiver-examples/main/notebooks/_static/street_elev.png
    :target: https://github.com/cheginit/HyRiver-examples/blob/main/notebooks/3dep.ipynb
    :align: center

Contributing
------------

Contributions are very welcomed. Please read
`CONTRIBUTING.rst <https://github.com/cheginit/pygeoogc/blob/master/CONTRIBUTING.rst>`__
file for instructions.
