"""Command-line interface for Py3DEP."""
from pathlib import Path
from typing import List, Union

import click
import geopandas as gpd
import pandas as pd

from . import py3dep
from .exceptions import MissingColumns, MissingCRS
from .py3dep import DEF_CRS, LAYERS


def get_target_df(
    tdf: Union[pd.DataFrame, gpd.GeoDataFrame], req_cols: List[str]
) -> Union[pd.DataFrame, gpd.GeoDataFrame]:
    """Check if all required columns exists in the dataframe.

    It also re-orders the columns based on ``req_cols`` order.
    """
    missing = [c for c in req_cols if c not in tdf]
    if len(missing) > 0:
        raise MissingColumns(missing)
    return tdf[req_cols]


CONTEXT_SETTINGS = {"help_option_names": ["-h", "--help"]}


@click.group(context_settings=CONTEXT_SETTINGS)
@click.pass_context
@click.option(
    "-q",
    "--query_source",
    default="tnm",
    type=click.Choice(["airmap", "tnm"], case_sensitive=False),
    help="Source of the elevation data.",
)
@click.option(
    "-s",
    "--save_dir",
    default="topo_3dep",
    type=click.Path(exists=False),
    help=(
        "Path to a directory to save the requested files."
        + "Extension for the outputs is either `.nc` for geometry or `.csv` for coords."
    ),
)
def cli(ctx: click.Context, query_source: str = "tnm", save_dir: Union[str, Path] = "topo_3dep"):
    """Command-line interface for Py3DEP."""
    Path(save_dir).mkdir(parents=True, exist_ok=True)
    ctx.obj = {"query_source": query_source, "save_dir": save_dir}


@cli.command("coords", context_settings=CONTEXT_SETTINGS)
@click.pass_context
@click.argument("fpath", type=click.Path(exists=True))
@click.argument("crs", type=str)
def coords(
    ctx: click.Context,
    fpath: Path,
    crs: str,
):
    r"""Retrieve topographic data for a list of coordinates.

    FPATH: Path to a csv file with two columns named ``x`` and ``y``.

    CRS: CRS of the input coordinates.

    Examples:

        $ py3dep  -s topo_dir coords ny_coords.csv epsg:4326
    """  # noqa: D412
    elev = get_target_df(pd.read_csv(fpath), ["x", "y"])

    item_str = "items" if len(elev) > 1 else "item"
    click.echo(f"Found {len(elev)} {item_str} in {fpath}. Retrieving ...")
    coords = list(elev.itertuples(index=False, name=None))
    elev["elevation"] = py3dep.elevation_bycoords(coords, crs, ctx.obj["query_source"])
    elev.astype("f8").to_csv(Path(ctx.obj["save_dir"], f"{Path(fpath).stem}_elevation.csv"))


@cli.command("geometry", context_settings=CONTEXT_SETTINGS)
@click.pass_context
@click.argument("fpath", type=click.Path(exists=True))
@click.argument("layer", type=click.Choice(LAYERS, case_sensitive=True))
def geometry(
    ctx: click.Context,
    fpath: Path,
    layer: str,
):
    r"""Retrieve topographic data within geometries.

    FPATH: Path to a geospatial file (any file that ``geopandas.read_file`` can open).

    This file should have three columns and contain ``crs`` attribute:

        - ``id``: Feature identifiers that py3dep uses as the output netcdf/csv filenames.

        - ``res``: Target resolution in meters.

        - ``geometry``: A Polygon or MultiPloygon.

    LAYER: A valid layer name when requesting for topographic data.

    Examples:

        $ py3dep -q airmap geometry ny_geom.gpkg "Slope Map"
    """  # noqa: D412
    target_df = gpd.read_file(fpath)
    if target_df.crs is None:
        raise MissingCRS
    crs = target_df.crs

    target_df = get_target_df(target_df, ["id", "res", "geometry"])
    args_list = (
        (g, r, Path(ctx.obj["save_dir"], f"{i}.nc"))
        for i, r, g in target_df.itertuples(index=False)
    )

    item_str = "items" if len(target_df) > 1 else "item"
    click.echo(f"Found {len(target_df)} {item_str} in {fpath}. Retrieving ...")
    with click.progressbar(
        args_list, label=f"Getting {layer} from 3DEP", length=len(target_df)
    ) as bar:
        for g, r, p in bar:
            py3dep.get_map(layer, g, r, geo_crs=crs, crs=DEF_CRS).to_netcdf(p)
