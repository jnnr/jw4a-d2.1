r"""
Use the prepared cutout of ERA5 model level data to prepare capacity factors for AWE,
both onshore and offshore.
"""
from operator import itemgetter

import atlite
import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr


def convert_wind_awe(ds, turbine):
    """
    Convert wind speeds for turbine to wind energy generation.
    """
    V, POW, hub_height, P = itemgetter("V", "POW", "hub_height", "P")(turbine)

    # wnd_hub = windm.extrapolate_wind_speed(ds, to_height=hub_height)
    wnd_hub = ds.wind_speed

    def _interpolate(da):
        return np.interp(da, V, POW / P)

    da = xr.apply_ufunc(
        _interpolate,
        wnd_hub,
        input_core_dims=[[]],
        output_core_dims=[[]],
        output_dtypes=[wnd_hub.dtype],
        dask="parallelized",
    )

    da.attrs["units"] = "MWh/MWp"
    da = da.rename("specific generation")

    return da


if __name__ == "__main__":
    if "snakemake" not in globals():
        from lib.helpers import mock_snakemake

        snakemake = mock_snakemake("prepare_capacity_factors")

    cutout = atlite.Cutout(snakemake.input.cutout)

    # prepare power curve
    powercurve = pd.read_csv(snakemake.input.powercurve, index_col=0)
    powercurve = dict(
        V=np.array(list(powercurve.index)),
        POW=np.array(list(powercurve["power"])),
        hub_height=100,
    )
    powercurve["P"] = np.max(powercurve["POW"])

    # prepare capacity matrix
    # No need to define area per grid cell and capacity per area as we ask for
    # capacity factors per unit in the end
    if (
        snakemake.input.get("shapes") is not None
        and snakemake.input.get("availability") is not None
    ):
        raise ValueError("Either shapes or availability must be defined, not both.")

    elif snakemake.input.get("shapes") is not None:
        capacity_matrix = None
        index = None
        shapes = gpd.read_file(snakemake.input.shapes).set_index("id")

    elif snakemake.input.get("availability") is not None:
        shapes = None
        availability = xr.load_dataarray(snakemake.input.availability)
        capacity_matrix = availability.stack(spatial=["y", "x"])
        index = capacity_matrix.id

    # prepare capacity factors
    capfac_wind_awe = cutout.convert_and_aggregate(
        convert_func=convert_wind_awe,
        turbine=powercurve,
        per_unit=True,
        matrix=capacity_matrix,
        shapes=shapes,
        index=index,
    )

    capfac_wind_awe.to_netcdf(str(snakemake.output))
