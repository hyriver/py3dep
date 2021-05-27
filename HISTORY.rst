=======
History
=======

0.11.0 (unreleased)
-------------------

New Features
~~~~~~~~~~~~
- Added command-line interface (:issue_3dep:`10`).
- All feature query functions use persistent caching that can significantly improve
  the performance.

Breaking Changes
~~~~~~~~~~~~~~~~
- Drop support for python 3.6 since many of the dependencies have done so, such as
  ``xarray`` and ``pandas``.
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
