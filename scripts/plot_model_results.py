import pandas as pd
import plotnine as pn
import numpy as np



def aggregate_locs(df):
    _df = df.copy()
    _df["locs"] = _df["locs"].replace({"_\\d+$": ""}, regex=True)
    _df = _df.groupby(["locs", "techs"]).sum().reset_index().set_index("locs")
    return _df


def filter_techs(df, filter_techs):
    if not isinstance(filter_techs, list):
        filter_techs = [filter_techs]
    return df.loc[df.techs.apply(lambda x: any([f in x for f in filter_techs]), 1),:]


def drop_nan_inf(series):
    return series.loc[~series.isna() & (~series.isin([np.inf, -np.inf]))]


def plot_results(data, var_name, var_unit, sort=True):

    if sort:
        order = data.set_index(["locs", "techs"]).loc[:, var_name].groupby("locs").sum().sort_values(ascending=True).index
        data = data.assign(order=pd.Categorical(data["locs"], categories=order, ordered=True))

    plot = (
        pn.ggplot(data)
        + pn.geom_col(pn.aes(x="order", y=var_name, fill="techs"))
        + pn.labs(x="Region", y=f"{var_name} ({var_unit})")
        + pn.scale_color_discrete(guide=False)
    )

    return plot


if __name__ == "__main__":
    
    FACTOR_10000_MW_to_GW = 100
    RENEWABLES = ["open_field_pv", "roof_mounted_pv", "wind", "awe"]

    df = pd.read_csv(snakemake.input[0])

    df.loc[:, ["energy_cap", "energy_cap_max"]] = df.loc[:, ["energy_cap", "energy_cap_max"]] * FACTOR_10000_MW_to_GW

    df["share"] = df.energy_cap / df.energy_cap_max
    df.loc[df["energy_cap_max"] == np.inf, "share"] = np.nan

    df = aggregate_locs(df).reset_index()

    df_re = filter_techs(df, RENEWABLES)

    plot_results(df_re, "energy_cap", "GW").save(snakemake.output[0], dpi=300, height=6, width=12, facecolor="w", transparent=False)
