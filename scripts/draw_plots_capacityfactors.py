import xarray as xr
import pandas as pd
import plotnine as pn
import matplotlib.pyplot as plt
import geopandas as gpd

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


def plot_annual_average(capacity_factors, labels, fig=None, ax=None):
    if fig is None or ax is None:
        fig, ax = plt.subplots(figsize=(10, 5))

    df = capacity_factors.groupby("time.year").mean("time").to_dataframe()["__xarray_dataarray_variable__"].unstack("year")
    df = df.assign(mean=df.mean(axis=1)).sort_values('mean').drop('mean', axis=1)
    ticks = labels.loc[df.index]
    df.index = df.index.map(str)

    cmap = plt.get_cmap('rainbow', len(df.columns))
    for n, col in enumerate(df.columns):
        ax.plot(df[col], label=col, linestyle="none", marker="o", alpha=0.5, color=cmap(n))
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.set_xticklabels(ticks, rotation=90)
    plt.tight_layout()
    plt.grid(alpha=0.3)
    plt.ylim(0, 1)

    return fig, ax


def format_tuple(tupl):
    relevant_entries = [item for item in list(tupl) if isinstance(item, str)]
    result = "/".join(relevant_entries)
    return result


if __name__ ==  "__main__":
    if "snakemake" not in globals():
        from lib.helpers import mock_snakemake

        snakemake = mock_snakemake("draw_plots_capacityfactors")

    capacity_factors = xr.load_dataset(snakemake.input.path_capacity_factors)
    boundaries = gpd.read_file(snakemake.input.path_boundaries)

    # prepare labels to properly name the regions
    labels = boundaries[["iso_sov1", "iso_sov2"]].apply(format_tuple, axis=1)

    # plot annual average
    fig, ax = plot_annual_average(capacity_factors, labels)
    plt.savefig(snakemake.output.path_plot_average, dpi=300, transparent=False)

    # split index in year and month-day-hour
    df_capacity_factors = capacity_factors.to_dataframe()["__xarray_dataarray_variable__"].unstack("dim_0")
    df_capacity_factors["year"] = df_capacity_factors.index.year
    df_capacity_factors["month-day-hour"] = df_capacity_factors.index.strftime('%m %d %H')
    df_capacity_factors.set_index(["year", "month-day-hour"], inplace=True)

    # give proper names to regions
    df_capacity_factors.columns = labels.loc[df_capacity_factors.columns]
    df_capacity_factors.columns.name = "region"
    
    # sort data
    df = df_capacity_factors.reset_index()
    df = pd.pivot_table(df, index=df["month-day-hour"], columns="year", values=df.columns)
    df = sort_timeseries(df)

    # melt data to plot
    data = df
    melted = df.copy()
    melted = pd.DataFrame(melted.stack(["year", "region"]), columns=["var_value"])
    melted.index.names = ["sorted_hours", "year", "region"]
    melted = melted.reset_index()
    
    # sort regions by average capacity factor
    order = df.groupby("region", axis=1).mean().mean().sort_values().index.values
    melted["region"] = pd.Categorical(melted["region"], ordered=True, categories=order)

    # plot
    (
        pn.ggplot(melted, pn.aes(x="sorted_hours", y="var_value", color="year"))
        + pn.geom_line(alpha=0.5)
        + pn.facet_wrap("region", nrow=4, ncol=8)
        + pn.theme_minimal()
        + pn.labs(x="Sorted hours", y="Capacity factor", title="Load duration curves")
        + pn.theme(axis_text_x=pn.element_text(rotation=60, hjust=1))
    ).save(snakemake.output.path_plot, dpi=300, height=5, width=10, facecolor="w", transparent=False)
