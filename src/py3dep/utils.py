"""Utilities for Py3DEP."""

from __future__ import annotations

import functools
import heapq
import warnings
from typing import TYPE_CHECKING, Any, Callable, Literal, TypeVar, Union, overload

import numpy as np
import xarray as xr
from numpy.typing import NDArray

from py3dep.exceptions import InputTypeError, InputValueError, NoOutletError

if TYPE_CHECKING:
    import pyproj

    CRSTYPE = Union[int, str, pyproj.CRS]


__all__ = ["deg2mpm", "fill_depressions"]


try:
    from numba import config as numba_config
    from numba import njit, prange

    numba_config.THREADING_LAYER = "workqueue"
    has_numba = True
except ImportError:
    has_numba = False
    prange = range

    T = TypeVar("T")
    Func = Callable[..., T]

    def njit(
        signature_or_function: str | Func[T], nogil: bool = True, parallel: bool = False
    ) -> Callable[[Func[T]], Func[T]]:
        def decorator_njit(func: Func[T]) -> Func[T]:
            @functools.wraps(func)
            def wrapper_decorator(*args: tuple[Any, ...], **kwargs: dict[str, Any]) -> T:
                return func(*args, **kwargs)

            return wrapper_decorator

        return decorator_njit


FloatArray = NDArray[np.float32]
BoolArray = NDArray[np.bool_]
IntArray = NDArray[np.uint32]
DataArray = TypeVar("DataArray", FloatArray, xr.DataArray)


@njit(
    "boolean[:, ::1](f4[:, ::1], boolean[:, ::1], boolean[:, ::1], optional(uint32[::1]), optional(f4))",
    nogil=True,
    parallel=True,
)
def _get_queued(
    elevtn: FloatArray,
    not_done: BoolArray,
    structure: BoolArray,
    idxs_pit: IntArray | None,
    elv_max: np.float32 | None,
) -> BoolArray:
    """Get queued cells for depression filling."""
    nrow, ncol = elevtn.shape
    if idxs_pit is not None:
        # initiate queue with outlet cells
        queued = np.full((nrow, ncol), False, dtype=np.bool_)
        for idx in idxs_pit:
            queued.flat[idx] = True
        return queued
    # initiate queue with edge cells
    s = np.where(structure.ravel())[0]
    queued = not_done.copy()
    for r in prange(0, nrow):
        for c in prange(0, ncol):
            if ~not_done[r, c] or r in (0, nrow - 1) or c in (0, ncol - 1):
                continue
            a0 = not_done[slice(r - 1, r + 2), slice(c - 1, c + 2)].ravel()
            if np.all(a0[s]):
                queued[r, c] = False
    if elv_max is not None:
        queued = np.logical_and(queued, elevtn <= elv_max)
        if not np.any(queued):
            raise NoOutletError
    return queued


@njit(
    "f4[:, ::1](f4[:, ::1], unicode_type, optional(uint32[::1]), f4, f4, optional(f4), i4)",
    nogil=True,
)
def _fill_depressions(
    elevtn: FloatArray,
    outlets: Literal["edge", "min"],
    idxs_pit: IntArray | None,
    nodata: np.float32,
    max_depth: np.float32,
    elv_max: np.float32 | None,
    connectivity: np.uint8,
) -> FloatArray:
    """Fill local depressions in elevation data based on Wang and Liu (2006)."""
    nrow, ncol = elevtn.shape
    delv = elevtn.copy()
    done = np.isnan(elevtn) if np.isnan(nodata) else np.isclose(elevtn, nodata)

    struct = np.full((3, 3), True, dtype=np.bool_)
    if connectivity == 4:
        struct[0, 0], struct[-1, -1] = False, False
        struct[0, -1], struct[-1, 0] = False, False

    queued = _get_queued(elevtn, ~done, struct, idxs_pit, elv_max)

    # queue contains (elevation, row, col)
    # boundary is included to favor non-boundary cells over boundary
    # cells with same elevation
    q = [(np.float32(0), np.uint32(0), np.uint32(0)) for _ in range(0)]
    for r, c in zip(*np.where(queued)):
        heapq.heappush(q, (np.float32(elevtn[r, c]), np.uint32(r), np.uint32(c)))
    heapq.heapify(q)

    # restrict queue to global edge minimum (single outlet)
    if outlets == "min":
        q = [heapq.heappop(q)]
        queued = np.full((nrow, ncol), False, dtype=np.bool_)
        queued[q[0][-2], q[0][-1]] = True

    # loop over cells and neighbors with ascending cell elevation.
    drs, dcs = np.where(struct)
    drs, dcs = drs - 1, dcs - 1
    while q:
        z0, r0, c0 = heapq.heappop(q)
        for dr, dc in zip(drs, dcs):
            r = r0 + dr
            c = c0 + dc
            if not (0 <= r < nrow and 0 <= c < ncol) or done[r, c]:
                continue
            z1 = elevtn[r, c]
            # local depression if dz > 0
            dz = z0 - z1
            # if positive max_depth: don't fill when dz > max_depth
            if max_depth >= 0 and dz >= max_depth:
                if not queued[r, c]:
                    heapq.heappush(q, (np.float32(z1), np.uint32(r), np.uint32(c)))
                    queued[r, c] = True
                continue
            if dz > 0:
                delv[r, c] += dz
            # add to queue if not already in queue
            if ~queued[r, c]:
                heapq.heappush(q, (np.float32(delv[r, c]), np.uint32(r), np.uint32(c)))
                queued[r, c] = True
            done[r, c] = True
    return delv


def fill_depressions(
    elevtn: DataArray,
    outlets: Literal["edge", "min"] = "min",
    idxs_pit: IntArray | None = None,
    nodata: float = np.nan,
    max_depth: float = -1.0,
    elv_max: float | None = None,
    connectivity: Literal[4, 8] = 8,
) -> DataArray:
    """Fill local depressions in elevation data based on Wang and Liu (2006).

    .. note::

        This function is based on the ``fill_depressions`` function from the
        `pyflwdir <https://github.com/Deltares/pyflwdir>`__ package. This function
        improves the performance of the original function by a factor of up to 2 and
        adds more input checks. Additionally, it works with ``xarray.DataArray`` objects.

        Outlets are assumed to occur at the edge of valid elevation cells ``outlets='edge'``;
        at the lowest valid edge cell to create one single outlet ``outlets='min'``;
        or at user provided outlet cells ``idxs_pit``.

        Depressions elsewhere are filled based on its lowest pour point elevation.
        If the pour point depth is larger than the maximum pour point depth ``max_depth``
        a pit is set at the depression local minimum elevation.

        Wang, L., & Liu, H. (2006). https://doi.org/10.1080/13658810500433453

    Parameters
    ----------
    elevtn: numpy.ndarray or xarray.DataArray
        elevation raster as a 2D ``numpy.ndarray`` or ``xarray.DataArray``.
    outlets: {"edge", "min}, optional
        Initial basin outlet(s) at the edge of all cells ('edge')
        or only the minimum elevation edge cell ('min'; default)
    idxs_pit: 1D array of int, optional
        Linear indices of outlet cells, in any, defaults to None.
    nodata: float, optional
        nodata value, defaults to ``numpy.nan``.
    max_depth: float, optional
        Maximum pour point depth. Depressions with a larger pour point
        depth are set as pit. A negative value (default) equals an infitely
        large pour point depth causing all depressions to be filled.
        Defaults to -1.0.
    elv_max, float, optional
        Maximum elevation for outlets, only in combination with ``outlets='edge'``.
        By default ``None``.
    connectivity: {4, 8}, optional
        Number of neighboring cells to consider, defaults to 8.

    Returns
    -------
    elevtn_out: numpy.ndarray
        Depression filled elevation with type float32.
    """
    if not has_numba:
        warnings.warn(
            "Numba not installed. Using very slow pure python version.", UserWarning, stacklevel=2
        )
    if not isinstance(elevtn, (np.ndarray, xr.DataArray)):
        raise InputTypeError("elevtn", "2D numpy.ndarray or xarray.DataArray")
    if elevtn.ndim != 2:
        raise InputTypeError("elevtn", "2D numpy.ndarray or xarray.DataArray")
    if outlets not in ("edge", "min"):
        raise InputValueError("outlets", ("edge", "min"))
    if connectivity not in (4, 8):
        raise InputValueError("connectivity", (4, 8))
    if idxs_pit is not None:  # noqa: SIM102
        if not isinstance(elevtn, np.ndarray) and idxs_pit.ndim != 1:
            raise InputTypeError("idxs_pit", "1D numpy.ndarray")
    _elevtn = np.asarray(elevtn, dtype=np.float32)
    _idxs_pit = None if idxs_pit is None else np.asarray(idxs_pit, dtype=np.uint32)
    _nodata = np.float32(nodata)
    _max_depth = np.float32(max_depth)
    _elv_max = None if elv_max is None else np.float32(elv_max)
    _connectivity = np.uint8(connectivity)
    corrected = _fill_depressions(
        _elevtn, outlets, _idxs_pit, _nodata, _max_depth, _elv_max, _connectivity
    )
    if isinstance(elevtn, xr.DataArray):
        return elevtn.copy(data=corrected)
    return corrected


def deg2mpm(slope: xr.DataArray) -> xr.DataArray:
    """Convert slope from degrees to meter/meter.

    Parameters
    ----------
    slope : xarray.DataArray
        Slope in degrees.

    Returns
    -------
    xarray.DataArray
        Slope in meter/meter. The name is set to ``slope`` and the ``units`` attribute
        is set to ``m/m``.
    """
    with xr.set_options(keep_attrs=True):
        nodata = slope.rio.nodata
        nodata = np.nan if nodata is None else nodata
        if not np.isnan(nodata):
            slope = slope.where(slope != nodata, drop=False)
        slope = xr.where(slope == 90, np.nan, slope)
        slope = np.tan(np.deg2rad(slope))  # pyright: ignore[reportAssignmentType]
        slope = slope.rio.write_nodata(np.nan)
        slope.name = "slope"
        slope.attrs["units"] = "m/m"
    return slope


@overload
def rename_layers(ds: xr.DataArray, valid_layers: list[str]) -> xr.DataArray: ...


@overload
def rename_layers(ds: xr.Dataset, valid_layers: list[str]) -> xr.Dataset: ...


def rename_layers(
    ds: xr.DataArray | xr.Dataset, valid_layers: list[str]
) -> xr.DataArray | xr.Dataset:
    """Rename layers in a dataset."""
    rename = {lyr: lyr.split(":")[-1].replace(" ", "_").lower() for lyr in valid_layers}
    rename.update({"3DEPElevation:None": "elevation"})

    if isinstance(ds, xr.DataArray):
        ds.name = rename[str(ds.name)]
    else:
        ds = ds.rename({n: rename[str(n)] for n in ds})
    return ds
