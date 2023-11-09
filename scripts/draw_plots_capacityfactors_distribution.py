import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr


def read_multiple_ncfiles(paths: dict) -> pd.DataFrame:
    capacity_factors = {}
    for name, path_capacity_factors in paths.items():
        ds = xr.load_dataset(path_capacity_factors)
        ds = ds.rename({"__xarray_dataarray_variable__": name})
        capacity_factors[name] = ds
    return capacity_factors


def plot_boxplot(df: pd.DataFrame, ax: plt.Axes = None) -> plt.Axes:
    if ax is None:
        fig, ax = plt.subplots()

    pass


def plot_histogram():
    pass


if __name__ == "__main__":
    snakemake_inputs = dict(
        offshore_deep_awe="build/capacity_factors/capacity_factors_offshore_deep_awe.nc",
        offshore_shallow_awe="build/capacity_factors/capacity_factors_offshore_shallow_awe.nc",
        old_wind_offshore="build/capacity_factors/capacity_factors_old_wind-offshore.nc",
        old_wind_onshore="build/capacity_factors/capacity_factors_old_wind-onshore.nc",
        onshore_awe="build/capacity_factors/capacity_factors_onshore_awe.nc",
    )
    snakemake_outputs = "build/plots/capacityfactors_distribution.png"

    # Read data
    capacity_factors = read_multiple_ncfiles(snakemake_inputs)
    import pdb

    pdb.set_trace()
    print(capacity_factors)
