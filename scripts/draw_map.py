r"""
Inputs
------
- water_depth
- boundaries
- potentials_onshore

Outputs
-------

Description
-----------
Plot maps of potentials

Plot potential cost curves
"""
import xarray as xr
import geopandas as gpd
import matplotlib.pyplot as plt
import rioxarray
from matplotlib.colors import ListedColormap


def map_area_potential(water_depth, boundaries_eez, boundaries_onshore, natura2000, fig=None, ax=None):
    r"""
    Shows exclusion criteria and regions for onshore and offshore wind, deep and shallow
    """
    if fig is None or ax is None:
        fig, ax = plt.subplots(figsize=(12, 10))

    water_depth = water_depth.coarsen(dim={"lat": 20, "lon": 20}, boundary="trim").mean()
    water_depth = water_depth.rio.write_crs("epsg:4326")
    water_depth = water_depth.rio.clip(boundaries_eez.geometry)
    water_depth = water_depth.drop("spatial_ref")

    # create masks for shallow and deep water
    threshold = -60
    water_shallow = ((water_depth.elevation < 0) & (water_depth.elevation > threshold)).astype(int)
    water_deep = ((water_depth.elevation < 0) & (water_depth.elevation < threshold)).astype(int)

    # define colormaps
    colors = [(0, 0, 0, 0), "#1aeaef"]
    cmap1 = ListedColormap(colors)
    colors = [(0, 0, 0, 0), "#675fb5"]
    cmap2 = ListedColormap(colors)

    # TODO: Crete and Kosovo are missing and some region in front of the coast of Normandy
    water_shallow.plot(ax=ax, cmap=cmap1, add_colorbar=False)
    water_deep.plot(ax=ax, cmap=cmap2, add_colorbar=False)
    boundaries_eez.geometry.boundary.plot(ax=ax, linewidth=0.5, color='#000000')
    boundaries_onshore.geometry.boundary.plot(ax=ax, linewidth=0.5, color='#100000')
    boundaries_onshore.geometry.plot(ax=ax, color='#000000', alpha=0.4)
    natura2000.geometry.to_crs("epsg:4326").plot(ax=ax, color='#369249', alpha=1)

    plt.tight_layout()
    ax.set_axis_off()

    return fig, ax


def map_wind_speeds(wind_speed, boundaries, fig=None, ax=None):
    r"""
    Shows capacityfactors for both onshore and offshore, can choose between hawt and awe
    """
    if fig is None or ax is None:
        fig, ax = plt.subplots(figsize=(12, 10))

    # prepare mean wind speed
    wind_speed_mean = wind_speed.wind_speed.mean("time")
    wind_speed_mean = wind_speed_mean.rio.write_crs("epsg:4326").drop("spatial_ref")

    # plot
    boundaries.geometry.boundary.plot(ax=ax, linewidth=0.5, color='#000000')
    wind_speed_mean.plot(cmap="Blues", ax=ax)

    ax.set_axis_off()

    return fig, ax


if __name__ == "__main__":
    if "snakemake" not in globals():
        from lib.helpers import mock_snakemake

        snakemake = mock_snakemake("draw_map")

    # load and prepare data
    natura2000 = gpd.read_file(snakemake.input.natura2000)

    boundaries_eez = gpd.read_file(snakemake.input.boundaries_eez)
    boundaries_onshore = gpd.read_file(snakemake.input.boundaries_onshore)

    water_depth = xr.open_dataset(snakemake.input.water_depth)

    wind_speed_era5_model_level = xr.load_dataset(snakemake.input.cutout)
    
    map_area_potential(water_depth, boundaries_eez, boundaries_onshore, natura2000)
    plt.savefig(snakemake.output.areas, dpi=300, bbox_inches="tight", transparent=False)

    map_wind_speeds(wind_speed_era5_model_level, boundaries_eez)
    plt.savefig(snakemake.output.wind_speeds, dpi=300, bbox_inches="tight", transparent=False)
