=======
History
=======

0.2.1 (unreleased)
------------------

- Add support for saving maps as ``geotiff`` file(s).

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

- Add a new function, ``show_versions``, for getting versions of the installed dependecies
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
