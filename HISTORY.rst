=======
History
=======

0.13.0 (2022-04-03)
-------------------

New Features
~~~~~~~~~~~~
- Add a new function called ``query_3dep_sources`` for querying bounds of 3DEP's
  data sources within a bounding box. It returns a geo-dataframe that contains
  the bounding box of each data source and a column ``dem_res`` identifying the
  resolution of the raw topographic data within each geometry.
- Add a new function called ``elevation_profile`` for getting elevation profile
  along a line at a given spacing. This function converts the line to a B-spline
  and then calculates the elevation along the spline at a given uniform spacing.

Breaking Changes
~~~~~~~~~~~~~~~~
- Remove caching-related arguments from all functions since now they
  can be set globally via three environmental variables:

  * ``HYRIVER_CACHE_NAME``: Path to the caching SQLite database.
  * ``HYRIVER_CACHE_EXPIRE``: Expiration time for cached requests in seconds.
  * ``HYRIVER_CACHE_DISABLE``: Disable reading/writing from/to the cache file.

  You can do this like so:

.. code-block:: python

    import os

    os.environ["HYRIVER_CACHE_NAME"] = "path/to/file.sqlite"
    os.environ["HYRIVER_CACHE_EXPIRE"] = "3600"
    os.environ["HYRIVER_CACHE_DISABLE"] = "true"

0.12.2 (2022-01-15)
-------------------

New Features
~~~~~~~~~~~~
- Add a new DEM source to ``elevation_bycoords`` to get elevation from
  the National Map's 3DEP WMS service. This can replace the ``tnm`` source
  since ``tnm`` is not stable.
- Add a new function called ``check_3dep_availability`` to check the availability
  of 3DEP's native resolutions within an area of interest. It returns a ``dict``
  with keys corresponding to the available resolutions and its values are boolean
  values indicating whether the resolution is available or not.
- Replace no data values of ``slope`` in ``deg2mm`` with ``np.nan``, so they do not
  get converted to another value. The output of this function has ``np.float64`` type.

Internal Changes
~~~~~~~~~~~~~~~~
- Refactor ``ElevationByCoords`` by using ``__post_init__`` for validating the
  input parameters rather than ``pydantic``'s validators.
- Refactor ``elevation_bygrid`` by using ``get_map`` to get DEM and ``rioxarray``
  for re-projection.
- Add type checking with ``typeguard`` and fixed typing issues raised by
  ``typeguard``.
- Refactor ``show_versions`` to ensure getting correct versions of all
  dependencies.

0.12.1 (2021-12-31)
-------------------

Internal Changes
~~~~~~~~~~~~~~~~
- Use the three new ``ar.retrieve_*`` functions instead of the old ``ar.retrieve``
  function to improve type hinting and to make the API more consistent.

0.12.0 (2021-12-27)
-------------------

Breaking Changes
~~~~~~~~~~~~~~~~
- Set the request caching's expiration time to never expire. Add two flags to all
  functions to control the caching: ``expire_after`` and ``disable_caching``.

Internal Changes
~~~~~~~~~~~~~~~~
- Add all the missing types so ``mypy --strict`` passes.
- Improve performance of ``elevation_bygrid`` by ignoring unnecessary validation.

0.11.4 (2021-11-12)
-------------------

Internal Changes
~~~~~~~~~~~~~~~~
- Use ``rioxarray`` for dealing with ``GeoTIFF`` binaries since ``xarray``
  deprecated the ``xarray.open_rasterio`` function, as it's discussed
  in this `PR <https://github.com/pydata/xarray/pull/5808>`__.
- Use ``importlib-metadata`` for getting the version instead of ``pkg_resources``
  to decrease import time as discussed in this
  `issue <https://github.com/pydata/xarray/issues/5676>`__.

0.11.3 (2021-10-03)
-------------------

Breaking Changes
~~~~~~~~~~~~~~~~
- Rewrite the command-line interface using ``click.group`` to improve UX.
  The command is now ``py3dep [command] [args] [options]``. The two supported commands are
  ``coords`` for getting elevations of a dataframe of coordinates in ``EPSG:4326`` CRS
  and ``geometry`` for getting the elevation of a geo-dataframe of geometries. Each sub-command
  now has a separate help message. The format of the input file for the ``coords`` command
  is now ``csv`` and for the ``geometry`` command is ``.shp`` or ``.gpkg`` and must have a
  ``crs`` attribute. Also, the ``geometry`` command now accepts multiple layers via the
  ``--layers`` (``-l``) option. More information and examples can be in the ``README.rst`` file.

New Features
~~~~~~~~~~~~
- Make ``fill_depressions`` function public. This function conditions an input DEM
  by applying
  `depression filling <https://richdem.readthedocs.io/en/latest/depression_filling.html>`__
  and
  `flat area resolution <https://richdem.readthedocs.io/en/latest/flat_resolution.html>`__
  operations.

Internal Changes
~~~~~~~~~~~~~~~~
- The ``get_map`` function now checks for validation of the input ``layers`` argument before
  sending the actual request with a more helpful message.
- Improve docstrings.
- Move ``deg2mpm``, ``fill_depressions``, and ``reproject_gtiff`` functions to a new file
  called ``utils``. Both ``deg2mpm`` and ``fill_depressions`` functions are still accessible
  from ``py3dep`` directly.
- Increase the test coverage.
- Use one of the ``click``'s internal functions, ``click..testing.CliRunner``,
  to run the CLI tests.

0.11.2 (2021-09-17)
-------------------

Bug Fixes
~~~~~~~~~
- Fix a bug related to ``elevation_bycoords`` where CRS validation fails if its
  type is ``pyrpoj.CRS`` by converting inputs with CRS types to string.

Internal Changes
~~~~~~~~~~~~~~~~
- Fix a couple of typing issues and update the ``get_transform`` API based on the
  recent changes in ``pygeoutils`` v0.11.5.


0.11.1 (2021-07-31)
-------------------

The first highlight of this release is a major refactor of ``elevation_bycoords`` by
adding support for the Bulk Point Query Service and improving the overall performance
of the function. Another highlight is support for performing depression filling
in ``elevation_bygrid`` before sampling the underlying DEM.

New Features
~~~~~~~~~~~~
- Refactor ``elevation_bycoords`` function to add support for getting
  elevations of a list of coordinates via The National Map's
  `Point Query Service <https://apps.nationalmap.gov/bulkpqs/>`__. This service is more
  accurate than Airmap, but it's limited to the US only. You can select the source via
  a new argument called ``source``. You can set it to ``source=tnm`` to use the TNM
  service. The default is ``tnm``.
- Refactor ``elevation_bygrid`` function to add a new capability via ``fill_depressions``
  argument for filling depressions in the obtained DEM before extracting elevation data
  for the input grid points. This is achieved via
  `RichDEM <https://richdem.readthedocs.io>`__ that needs to be installed if this
  functionality is desired. You can install it via ``pip`` or ``conda`` (``mamba``).

Internal Changes
~~~~~~~~~~~~~~~~
- Migrate to using ``AsyncRetriever`` for handling communications with web services.
- Handle the interpolation step in ``elevation_bygrid`` function more efficiently
  using ``xarray``.

0.11.0 (2021-06-19)
-------------------

New Features
~~~~~~~~~~~~
- Added command-line interface (:issue_3dep:`10`).
- All feature query functions use persistent caching that can significantly improve
  the performance.

Breaking Changes
~~~~~~~~~~~~~~~~
- Drop support for Python 3.6 since many of the dependencies such as ``xarray`` and ``pandas``
  have done so.
- The returned ``xarray`` objects are in parallel mode, i.e., in some cases ``compute`` method
  should be used to get the results.
- Save the output as a ``netcdf`` instead of ``raster`` since conversion
  from ``nc`` to ``tiff`` can be easily done with ``rioxarray``.

0.10.1 (2021-03-27)
-------------------

- Add announcement regarding the new name for the software stack, HyRiver.
- Improve ``pip`` installation and release workflow.

0.10.0 (2021-03-06)
-------------------

- The first release after renaming hydrodata to PyGeoHydro.
- Make ``mypy`` checks more strict and fix all the errors and prevent possible
  bugs.
- Speed up CI testing by using ``mamba`` and caching.

0.9.0 (2021-02-14)
------------------

- Bump version to the same version as PyGeoHydro.
- Add support for saving maps as ``geotiff`` file(s).
- Replace ``Elevation Point Query Service`` service with ``AirMap`` for getting
  elevations for a list of coordinates in bulk since ``AirMap`` is much faster.
  The resolution of ``AirMap`` is 30 m.
- Use ``cytoolz`` for some operations for improving performance.

0.2.0 (2020-12-06)
------------------

- Add support for multipolygon.
- Remove the ``fill_hole`` argument.
- Add a new function to get elevations for a list of coordinates called ``elevation_bycoords``.
- Refactor ``elevation_bygrid`` function for increasing readability and performance.

0.1.7 (2020-08-18)
------------------

- Added a rename operation to ``get_map`` to automatically rename the variables to a
  more sensible one.
- Replaced ``simplejson`` with ``orjson`` to speed-up JSON operations.

0.1.6 (2020-08-11)
------------------

- Add a new function, ``show_versions``, for getting versions of the installed dependencies
  which is useful for debugging and reporting.
- Fix typos in the docs and improved the README.
- Improve testing and coverage.

0.1.5 (2020-08-03)
------------------

- Fixed the geometry CRS issue
- Improved the documentation

0.1.4 (2020-07-23)
------------------

- Refactor ``get_map`` to use ``pygeoutils`` package.
- Change the versioning method to ``setuptools_scm``.
- Polish README and add installation from ``conda-forge``.

0.1.0 (2020-07-19)
------------------

- First release on PyPI.
