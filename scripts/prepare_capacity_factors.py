from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr


def read_power_curve(filepath_powercurve: Path or str) -> xr.Dataset:
    r"""
    Read power curve from csv file.
    """
    power_curve = pd.read_csv(filepath_powercurve, index_col="windspeed")
    power_curve = xr.Dataset.from_dataframe(power_curve)

    return power_curve


def join_netcdfs(datasets: list[xr.Dataset]) -> xr.Dataset:
    r"""
    Joins a list of xarray datasets into a single dataset.
    """
    dataset_joined = xr.concat(datasets, dim="location")

    return dataset_joined


def read_windspeed(filepaths_windspeed: list[Path or str]) -> xr.Dataset:
    r"""
    Read wind speed data from several netcdf files.
    """
    list_windspeed = []
    for path in filepaths_windspeed:
        data = xr.open_dataset(path)
        data = data.assign_coords(location=path.stem)
        data = (
            data.reset_index(["latitude", "longitude"])
            .reset_coords(["latitude", "longitude"])
            .squeeze()
        )
        list_windspeed.append(data)

    windspeed = join_netcdfs(list_windspeed)

    return windspeed


def calculate_capacity_factors(
    windspeed: xr.Dataset, power_curve: xr.Dataset
) -> xr.Dataset:
    r"""
    Calculate capacity factors from wind speed data and a power curve.

    Parameters
    ----------
    windspeeds : :class:`xarray.Dataset`
        Wind speed data with a datetime index.
    powercurve : :class:`xarray.Dataset`

    Returns
    -------
    capacity_factors : :class:`xarray.Dataset`
        Capacity factors with a datetime index.
    """
    data = windspeed.copy()

    # calculate absolute wind speed
    data = data.assign(windspeed=lambda x: np.sqrt(x.u**2 + x.v**2))  # pythagoras

    # round to 1 decimal to be able to merge with power curve
    data = data.assign(windspeed=lambda x: np.round(x.windspeed, decimals=1))

    # apply power curve
    capacity_factors = power_curve["capfac"].sel(
        windspeed=data.windspeed, method="nearest"
    )

    return capacity_factors


if __name__ == "__main__":
    if "snakemake" not in globals():
        import sys
        from pathlib import Path

        sys.path.append(str(Path(__file__).parent.parent))
        from lib.helpers import mock_snakemake

        snakemake = mock_snakemake("prepare_capacity_factors")

    filepaths_windspeed = Path(snakemake.input.windspeed).glob("*.nc")
    filepath_powercurve = snakemake.input.powercurve
    filepath_capacity_factors = snakemake.output.capacity_factors

    power_curve = read_power_curve(filepath_powercurve)
    windspeed = read_windspeed(filepaths_windspeed)

    capacity_factors = calculate_capacity_factors(windspeed, power_curve)

    df_capacity_factors = (
        capacity_factors.to_dataframe()
        .reset_index()
        .pivot(index="time", columns="location", values="capfac")
    )

    df_capacity_factors.to_csv(filepath_capacity_factors)
