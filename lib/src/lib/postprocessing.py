import calliope
import pandas as pd
from pathlib import Path


here = Path(__file__).parent
CONVERSION_FACTORS = pd.read_csv(here / "unit_conversion.csv", index_col=0).loc["converion_factor"].to_dict()


def get_formatted_data(model: calliope.Model, coefficient: str, conversion_factors=None) -> pd.Series:
    r"""
    Get input or results data from a calliope model.
    
    Parameters
    ----------
    model : calliope.Model
        A calliope model.
    coefficient : str
        The coefficient to get from the model.
    conversion_factors : dict (default: None)
        Dictionary mapping the coefficient to a conversion factor.
    
    Returns
    -------
    data : pd.Series
        Requested data.
    """
    if conversion_factors is None:
        conversion_factors = CONVERSION_FACTORS


    data = (
        model.get_formatted_array(coefficient, index_format="multiindex")
        .to_dataframe()
    )

    data *= conversion_factors[coefficient]

    return data


def aggregate(data: pd.Series or pd.DataFrame, agg_map: dict, axis=0: int) -> pd.Series or pd.DataFrame:
    r"""
    Aggregate data 
    
    Parameters
    ----------
    data : pd.Series or pd.DataFrame
        Data to aggregate.
    agg_map : dict
        Dictionary mapping the indices/columns to aggregate.
    axis : int (default: 0)
        Axis to aggregate along.

    Returns
    -------
    data_agg : pd.Series or pd.DataFrame
        Aggregated data.
    """
    data_agg = data.copy()

    data_agg = data_agg.groupby(["locs", "techs"]).sum().reset_index().set_index("locs")

    return data_agg


def filter_index(data: pd.Series or pd.DataFrame, where: str or list[str], level: str) -> pd.Series or pd.DataFrame:
    r"""
    Filter DataFrame by index.
    """
    pass


def df_drop_nan_inf(data: pd.DataFrame, axis: int, drop_if_any=True: bool) -> pd.DataFrame:
    r"""
    Drop nan and inf values from a dataframe.

    Parameters
    ----------
    data : pd.DataFrame
        Dataframe to clean.
    axis : int
        Axis to clean along.
    drop_if_any : bool (default: True)
        If True, drop if any value is nan or inf.
        If False, drop if all values are nan or inf.

    Returns
    -------
    pd.DataFrame
        Cleaned dataframe.
    """
    if drop_if_any:
        condition = ~(df.isna() | df.isin([np.inf, -np.inf])).any(axis)
    else:
        condition = ~(df.isna() | df.isin([np.inf, -np.inf])).all(axis)
    if axis == 1:
        return df.loc[condition, :]
    elif axis == 0:
        return df.loc[:, condition]
    else:
        raise ValueError("axis must be 0 or 1")
