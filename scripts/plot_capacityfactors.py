import xarray as xr
import pandas as pd
import plotnine as pn


def melt_data(data):
    r"""
    Melt data to be able to plot it with plotnine.

    Input
    -----
    data : pandas.DataFrame
        Data to be melted.

    Returns
    -------
    melted : pandas.DataFrame
        Melted data.

    Examples
    --------

    Melts a DataFrame with a timeindex and multiple columns

        timeindex  |  column_1  |  ...
        ----------------------------
        2023-08-31 |    0.1     |  ...
        ...

    to a DataFrame with index 'id' and columns 'timeindex', 'var_name' and 'var_value':

        id | timeindex  | var_name | var_value
        ----------------------------
        0  | 2023-08-31 | column_1  | 0.1
        ...
    """
    melted = data.reset_index()
    melted = pd.melt(melted, id_vars=data.index.name, value_vars=data.columns, var_name="var_name", value_name="var_value")
    melted.index.name = "id"

    return melted


def sort_timeseries(data, ascending=True):
    r"""
    Sorts each of the columns of a DataFrame separately.
    Drops index and resets it to a range from 0 to len(data).

        index |  column_1  |  ...
        ----------------------------
        0     |  0.1       |   ...
        ...
    """
    sorted_data = data.reset_index(drop=True)
    sorted_data.index.name = "id"

    for column in data.columns:
        sorted_data[column] = sorted_data[column].sort_values(ascending=ascending).values

    return sorted_data


if __name__ ==  "__main__":
    path_capacity_factors = "build/capacity_factors/capacity_factors_offshore_deep_awe.nc"
    path_plot = "build/plots/load_duration_wind_onshore_awe.png"

    capacity_factors = xr.load_dataset(path_capacity_factors)

    # split index in year and month-day-hour
    df_capacity_factors = capacity_factors.to_dataframe()["__xarray_dataarray_variable__"].unstack("dim_0")
    df_capacity_factors["year"] = df_capacity_factors.index.year
    df_capacity_factors["month-day-hour"] = df_capacity_factors.index.strftime('%m %d %H')
    df_capacity_factors.set_index(["year", "month-day-hour"], inplace=True)

    # sort data
    df = df_capacity_factors.reset_index()
    df = pd.pivot_table(df, index=df["month-day-hour"], columns="year", values=df.columns)
    df = sort_timeseries(df)

    # melt data to plot
    data = df
    melted = df.copy()
    melted = pd.DataFrame(melted.stack(["year", "dim_0"]), columns=["var_value"])
    # melted = pd.melt(melted, id_vars=data.index.names, value_vars=data.columns.to_list(), var_name="var_name", value_name="var_value")
    # melted.set_index(["year", "day-time", "var_name"], inplace=True)
    melted.index.names = ["sorted_hours", "year", "region"]
    melted = melted.reset_index()

    (
        pn.ggplot(melted, pn.aes(x="sorted_hours", y="var_value", color="region", group="region"))
        + pn.geom_line(alpha=0.5)
        + pn.facet_wrap("region", nrow=4, ncol=8)
        + pn.theme_minimal()
        + pn.labs(x="Sorted hours", y="Capacity factor", title="Load duration curves")
    ).save(path_plot, dpi=300, height=5, width=10, transparent=False)
