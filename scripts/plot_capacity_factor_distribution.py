import pandas as pd
import plotnine as pn
import xarray as xr


def read_multiple_ncfiles(paths: dict) -> pd.DataFrame:
    capacity_factors = {}
    for name, path_capacity_factors in paths.items():
        ds = xr.load_dataset(path_capacity_factors)
        ds = ds.rename({"__xarray_dataarray_variable__": name})
        capacity_factors[name] = ds.to_dataframe()

    capacity_factors = pd.concat(capacity_factors.values(), axis=1)

    return capacity_factors


def plot_boxplot(df: pd.DataFrame, facets: str = None) -> pn.ggplot:
    """
    Plot boxplot given a tidy dataframe.

    Parameters
    ----------
    df : pd.DataFrame
        Data to plot
    facets : str
        Variable to facet by (optional)

    Returns
    -------
    pn.ggplot
        The boxplot
    """
    boxplot = (
        pn.ggplot(df, pn.aes(x="variable", y="value"))
        + pn.geom_boxplot()
        + pn.facet_wrap(facets=facets)
        + pn.labs(x="Technology", y="Annual capacity factor")
        + pn.theme(axis_text_x=pn.element_text(angle=45))
    )

    return boxplot


def plot_histogram(
    df: pd.DataFrame, x: str, facets: str, nrow: int, ncol: int, bins: int
) -> pn.ggplot:
    """
    Plot histogram given a tidy dataframe.

    Parameters
    ----------
    df : pd.DataFrame
        Data to plot
    x : str
        _description_
    facets : str
        Variable to facet by
    nrow : int
        _description_
    ncol : int
        _description_
    bins : int
        _description_

    Returns
    -------
    pn.ggplot
        _description_
    """
    histogram = (
        pn.ggplot(df, pn.aes(x=x))
        + pn.facet_wrap(facets=facets)
        + pn.geom_histogram(alpha=0.5, bins=bins)
        + pn.labs(x="Annual capacity factor", y="Count")
        + pn.theme(legend_position="top")
    )

    return histogram


def df_drop_all_small(
    df: pd.DataFrame,
    eps: float,
    axis: int,
) -> pd.DataFrame:
    r"""
    Drop indices/columns with small values from a dataframe.

    Parameters
    ----------
    df : pd.DataFrame
        Dataframe to clean.
    eps : float
        Threshold for small values.
    axis : int
        Axis to clean along.

    Returns
    -------
    pd.DataFrame
        Cleaned dataframe.
    """
    condition = ~(abs(df) < eps).all(axis)
    if axis == 1:
        return df.loc[condition, :]
    elif axis == 0:
        return df.loc[:, condition]
    else:
        raise ValueError("axis must be 0 or 1")


if __name__ == "__main__":
    if "snakemake" not in globals():
        from lib.helpers import mock_snakemake

        snakemake = mock_snakemake("plot_capacity_factor_distribution")

    # Read data
    capacity_factors = read_multiple_ncfiles(snakemake.input)

    # Prepare data
    def drop_all_small_timeseries(df):
        df = df.unstack("id")
        df = df_drop_all_small(df, axis=0, eps=1e-3)
        df = df.stack("id")

        return df

    capacity_factors = drop_all_small_timeseries(capacity_factors)

    capacity_factors = capacity_factors.reset_index().melt(id_vars=["time", "id"])

    # Plot
    histogram = plot_histogram(
        capacity_factors, x="value", facet="variable", nrow=1, ncol=5, bins=30
    )
    histogram.save(
        snakemake.output.histogram,
        dpi=300,
        height=5,
        width=10,
        facecolor="w",
        transparent=False,
    )

    histogram = plot_histogram(
        capacity_factors,
        x="value",
        facets=("id", "variable"),
        nrow=len(capacity_factors.id.unique()),
        ncol=len(capacity_factors.variable.unique()),
        bins=30,
    )
    histogram.save(
        snakemake.output.regional_histogram,
        dpi=300,
        height=20,
        width=20,
        facecolor="w",
        transparent=False,
    )

    boxplot = plot_boxplot(capacity_factors, facets="id")
    boxplot.save(
        snakemake.output.regional_boxplot,
        dpi=300,
        height=20,
        width=20,
        facecolor="w",
        transparent=False,
    )
