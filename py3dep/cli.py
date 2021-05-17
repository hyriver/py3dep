"""Command-line interface for Py3DEP."""
import os
from pathlib import Path
from typing import List, Optional, Tuple, Union

import click
import geopandas as gpd
import pandas as pd
from shapely.geometry import MultiPolygon, Polygon

from . import py3dep
from .exceptions import MissingColumns


def get_target_df(
    tdf: Union[pd.DataFrame, gpd.GeoDataFrame], req_cols: List[str]
) -> Union[pd.DataFrame, gpd.GeoDataFrame]:
    """Check if all required columns exists in the dataframe.

    It also re-orders the columns based on req_cols order.
    """
    missing = [c for c in req_cols if c not in tdf]
    if len(missing) > 0:
        raise MissingColumns(missing)
    return tdf[req_cols]


def from_geometry(
    layer: str,
    geometry: Union[Polygon, MultiPolygon, Tuple[float, float, float, float]],
    res: float,
    crs: str,
    nc_path: Union[str, Path],
) -> None:
    """Get topographic data from 3DEP for a geometry."""
    py3dep.get_map(layer, geometry, res, geo_crs=crs, crs=crs, nc_path=nc_path)


def from_coords(coords: List[Tuple[float, float]], crs: str, csv_path: Union[str, Path]) -> None:
    """Get elevations of a set of coordinates in meter from airmap."""
    elev = pd.DataFrame.from_records(coords, columns=["x", "y"])
    elev["elevation"] = py3dep.elevation_bycoords(coords, crs)
    elev.astype("f8").to_csv(csv_path)


LAYERS = [
    "DEM",
    "Hillshade Gray",
    "Aspect Degrees",
    "Aspect Map",
    "GreyHillshade_elevationFill",
    "Hillshade Multidirectional",
    "Slope Map",
    "Slope Degrees",
    "Hillshade Elevation Tinted",
    "Height Ellipsoidal",
    "Contour 25",
    "Contour Smoothed 25",
]


@click.command()
@click.argument("target", type=click.Path(exists=True))
@click.argument("target_type", type=click.Choice(["geometry", "coords"], case_sensitive=False))
@click.argument("crs", type=str)
@click.option(
    "-l",
    "--layer",
    default=None,
    type=click.Choice(LAYERS, case_sensitive=True),
    help="Layer name when requesting for topographic data.",
)
@click.option(
    "-s",
    "--save_dir",
    type=click.Path(exists=False),
    default="topo_3dep",
    help="Path to a directory to save the requested files. Extension for the outputs is .nc for geometry and .csv for coords.",
)
def main(
    target: Path,
    target_type: str,
    crs: str,
    layer: Optional[str] = None,
    save_dir: Union[str, Path] = "topo_3dep",
):
    r"""Retrieve topographic data within geometries or elevations for a list of coordinates.

    TARGET: Path to a geospatial file (any file that geopandas.read_file can open) or a csv file.

    The geospatial file should have three columns:
    \n\t- id: Feature identifiers that py3dep uses as the output netcdf filenames.
    \n\t- res: Target resolution in meters.
    \n\t- geometry: A Polygon or MultiPloygon.

    The csv file should have two column: x and y.

    TARGET_TYPE: Type of input file: "coords" for csv, and "geometry" for geospatial.

    CRS: CRS of the input data.

    Example:

    \tpy3dep ny_coords.csv coords epsg:4326
    \n\tpy3dep ny_geom.gpkg geometry epsg:3857 --layer "Slope Map"
    """  # noqa: D412
    save_dir = Path(save_dir)
    target = Path(target)
    if not save_dir.exists():
        os.makedirs(save_dir, exist_ok=True)

    if target_type == "geometry":
        if layer is None:
            raise ValueError("layer option is required when target_type is geometry.")

        target_df = gpd.read_file(target, crs=crs)
        target_df = get_target_df(target_df, ["id", "res", "geometry"])
        args_list = (
            (layer, g, r, crs, Path(save_dir, f"{i}.nc"))
            for i, r, g in target_df.itertuples(index=False)
        )

        click.echo(f"Found {len(target_df)} items in {target}. Retrieving ...")
        with click.progressbar(
            args_list, label=f"Getting {layer} from 3DEP", length=len(target_df)
        ) as bar:
            for args in bar:
                from_geometry(*args)

        click.echo(f"Retrieved topography data for {len(target_df)} item(s).")
    else:
        target_df = pd.read_csv(target)
        target_df = get_target_df(target_df, ["x", "y"])
        click.echo(f"Found {len(target_df)} items in {target}. Retrieving ...")
        from_coords(
            list(target_df.to_records(index=False)),
            crs,
            Path(save_dir, f"{target.stem}_elevation.csv"),
        )
        click.echo(f"Retrieved elevation data for {len(target_df)} item(s).")
