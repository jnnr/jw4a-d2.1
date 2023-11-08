import pandas as pd
import plotnine as pn
import numpy as np


def stacked_bar_plot(df, x, y, fill, order=None, normalized=False, colors=None):
    r"""
    Plot a stacked bar plot.
    """
    _df = df.copy()

    if normalized:
        _df_max = _df.groupby(x).agg({y: "sum"}).reset_index()
        _df = _df.merge(_df_max, on=x, suffixes=("", "_max"))
        _df[y] = _df[y] / _df[y+"_max"]
        _df = _df.drop(columns=[y+"_max"])
        _df = _df.loc[~_df[y].isin([np.inf, -np.inf])]
    
    if colors:
        assert isinstance(colors, dict), "Colors should be a dictionary."
        assert _df.loc[:, fill].isin(colors.keys()).all(), "Not all colors are defined."
        color_map = pn.scale_fill_manual(breaks=list(colors.keys()), values=list(colors.values()))

    if order is not None:
        assert set(_df[x].unique()) == set(order), "Order should contain all items in x."
        categories = pd.Categorical(_df[x], categories=order, ordered=True)
        _df = _df.assign(**{x: categories})
    
    plot = (
        pn.ggplot(_df)

        + pn.geom_col(pn.aes(x=x, y=y, fill=fill))
        + pn.labs(x="Region", y=f"{var_name} ({var_unit})")
        + pn.scale_color_discrete(guide=False)
        + pn.theme(axis_text_x=pn.element_text(angle=90), legend_position="bottom")
        + color_map 
    )

    return plot