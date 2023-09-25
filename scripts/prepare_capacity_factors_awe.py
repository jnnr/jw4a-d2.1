r"""
Use the prepared cutout of ERA5 model level data to prepare capacity factors for AWE, both onshore and
offshore.
"""
import atlite
import functools
import numpy as np
import geopandas as gpd
import pandas as pd
import xarray as xr
from operator import itemgetter


def convert_wind_awe(ds, turbine):
    """
    Convert wind speeds for turbine to wind energy generation.
    """
    V, POW, hub_height, P = itemgetter("V", "POW", "hub_height", "P")(turbine)

    #wnd_hub = windm.extrapolate_wind_speed(ds, to_height=hub_height)
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
    # TODO https://atlite.readthedocs.io/en/latest/examples/landuse-availability.html

    boundaries = "build/shapes/eez.geojson"
    natura2000 = "data/potentials_offshore/natura2000_areas/eea_v_3035_100_k_natura2000_p_2021_v12_r01/SHP/Natura2000_end2021_rev1_epsg3035.shp"
    shipdensity = "data/potentials_offshore/shipping_density_global/shipdensity_global.tif"
    gebco = "data/potentials_offshore/gebco_2023_sub_ice_topo/GEBCO_2023_sub_ice_topo.nc"
    path_cutout = "build/cutouts/cutout-era5-model-level_adapted_copy.nc"
    path_powercurve = "data/power_curves/AWE_500kw_softwing.csv"
    snakemake_output = "build/capacity_factors/capacity_factors_awe.nc"

    df_boundaries = gpd.read_file(boundaries)
    cutout = atlite.Cutout(path_cutout)

    # prepare power curve
    awe_500kw_sw = pd.read_csv(path_powercurve, index_col=0)
    awe_500kw_sw = dict(V=np.array(list(awe_500kw_sw.index)), POW=np.array(list(awe_500kw_sw["power"])), hub_height=100)
    awe_500kw_sw["P"] = np.max(awe_500kw_sw["POW"])

    # prepare excluder
    excluder = atlite.ExclusionContainer()

    excluder.add_geometry(natura2000)

    max_depth = 60
    func = functools.partial(np.greater, -max_depth)
    excluder.add_raster(gebco, crs=4326, codes=func)

    # prepare availability matrix
    # availability = cutout.availabilitymatrix(df_boundaries, excluder)

    availability = xr.load_dataarray("build/capacity_factors/availability_awe.nc")
    # availability.to_netcdf("build/capacity_factors/availability_awe.nc")

    # prepare capacity matrix
    # No need to define area per grid cell and capacity per area as we ask for capacity factors per unit in the end
    capacity_matrix = availability.stack(spatial=["y", "x"])

    # prepare capacity factors
    capfac_wind_awe = cutout.convert_and_aggregate(convert_func=convert_wind_awe, turbine=awe_500kw_sw, per_unit=True, matrix=capacity_matrix)

    capfac_wind_awe.to_netcdf(snakemake_output)
