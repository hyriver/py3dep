=======
History
=======

0.11.3 (Unreleased)
-------------------

New Features
~~~~~~~~~~~~
- Rewrite the command line using ``click.group`` (command/sub-command) to improve UX.
  The command is now ``py3dep [command] [args] [options]``. Options are ``--save_dir``,
  (or ``-s``) and ``--query_source`` (or ``-q``). The two supported commands are
  ``coords`` for getting elevations of a list of coordinates and ``geometry`` for
  getting the elevation of within a geometry. Each sub-command now has a separate
  help message.

Internal Changes
~~~~~~~~~~~~~~~~
- The ``get_map`` function now checks for validation of the input ``layers`` argument before
  sending the actual request with a more helpful message.

0.11.2 (2021-09-17)
-------------------

Bug Fixes
~~~~~~~~~
- Fix a bug related to ``elevation_bycoords`` where crs validation fails if its
  type is ``pyrpoj.CRS`` by converting inputs with CRS types to string.

Internal Changes
~~~~~~~~~~~~~~~~
- Fix a couple of typing issues and update the ``get_transform`` API based on the
  recent changes in ``pygeoutils`` v0.11.5.


0.11.1 (2021-07-31)
-------------------

The first highlight of this release is a major refactor of ``elevation_bycoords`` by
adding support for the Bulk Point Query Service and improving overall performance
of the function. Another highlight is support for performing depression filling
in ``elevation_bygrid`` prior to sampling the underlying DEM.

New Features
~~~~~~~~~~~~
- Refactor ``elevation_bycoords`` function to add support for getting
  elevations of a list of coordinates via The National Map's
  `Point Query Service <https://apps.nationalmap.gov/bulkpqs/>`__. This service is more
  accurate than Airmap but it's limited to the US only. You can select the source via
  a new argument called ``source``. You can set it to ``source=tnm`` to use the TNM
  service. The default is ``tnm``.
- Refactor ``elevation_bygrid`` function to add a new capability via ``fill_depressions``
  argument for filling depressions in the obtained DEM before extracting elevation data
  for the input grid points.. This is achieved via
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

- Add announcement regarding the new name for the softwate stack, HyRiver.
- Improve ``pip`` installation and release workflow.

0.10.0 (2021-03-06)
-------------------

- The first release after renaming hydrodata to pygeohydro.
- Make ``mypy`` checks more strict and fix all the errors and prevent possible
  bugs.
- Speed up CI testing by using ``mamba`` and caching.

0.9.0 (2021-02-14)
------------------

- Bump version to the same version as pygeohydro.
- Add support for saving maps as ``geotiff`` file(s).
- Replace ``Elevation Point Query Service`` service with ``AirMap`` for getting
  elevations for a list of coordinates in bulk since ``AirMap`` is much faster.
  The resolution of ``AirMap`` is 30 m.
- Use ``cytoolz`` for some of the operations for improving performance.

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
