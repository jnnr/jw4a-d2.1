import calliope
import numpy as np
import pandas as pd

CONVERT_100_GW_to_GW = 100
CONVERT_10000km2_to_km2 = 10000
RENEWABLES = ["open_field_pv", "roof_mounted_pv", "wind", "awe"]


def get_tidy_data(model: calliope.Model, coefficient: str) -> pd.DataFrame:
    r"""
    Get data from a model and format it for plotting.
    """
    data = (
        model.get_formatted_array(coefficient, index_format="multiindex")
        .to_dataframe()
        .reset_index()
    )

    return data


def aggregate_locs(df):
    _df = df.copy()
    _df["locs"] = _df["locs"].replace({"_\\d+$": ""}, regex=True)
    _df = _df.groupby(["locs", "techs"]).sum().reset_index().set_index("locs")
    return _df


def filter_techs(df, filter_techs):
    if not isinstance(filter_techs, list):
        filter_techs = [filter_techs]
    return df.loc[df.techs.apply(lambda x: any([f in x for f in filter_techs]), 1), :]


def series_drop_nan_inf(series):
    return series.loc[~series.isna() & (~series.isin([np.inf, -np.inf]))]


def df_drop_nan_inf(df, axis, drop_if_any=True):
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


def prepare_overview(model, coefficient, techs, aggregate=False):
    df = get_tidy_data(model, coefficient)
    df = filter_techs(df, techs)
    if aggregate:
        df = aggregate_locs(df)
    df = df.set_index("techs", append=True).unstack("techs")
    df = df.reset_index()

    return df


def get_energy_cap(model):
    df = prepare_overview(model, "energy_cap_max", RENEWABLES, aggregate=True).round(2)
    df = df_drop_nan_inf(df, axis=0, drop_if_any=False).drop(
        "wind_offshore", axis=1, level=1
    )
    df *= CONVERT_100_GW_to_GW

    return df


def get_available_area_onshore(model):
    df = get_tidy_data(model, "available_area").round(2).set_index("locs")
    df = df_drop_nan_inf(df, axis=0, drop_if_any=False)
    df *= CONVERT_10000km2_to_km2

    return df


def get_offshore(model):
    selector = list(model.inputs.group_names_energy_cap_max.values)
    model.inputs.sel(group_names_energy_cap_max=selector)[
        ["group_energy_cap_max"]
    ].to_dataframe().unstack(1)
    df = get_tidy_data(model, "group_energy_cap_max").round(2)
    df["group"] = df["group_names_energy_cap_max"].str.extract(
        r"(wind_offshore_\w+_cap_max)"
    )
    df["locs"] = df["group_names_energy_cap_max"].str.extract(
        r"wind_offshore_\w+_cap_max_(.*)"
    )
    df = df.drop(columns=["group_names_energy_cap_max"])
    df = df.set_index(["locs", "group"]).unstack("group")

    return df


if __name__ == "__main__":
    if "snakemake" not in globals():
        from lib.helpers import mock_snakemake

        snakemake = mock_snakemake("table_area_constraints")

    path_inputs = snakemake.input[0]
    model = calliope.read_netcdf(path_inputs)

    # [x] wind_onshore_monopoly, rooftop_pv: energy_cap_max
    #       [ ] power_densities
    # [x] open_field_pv, wind_onshore_competing, awe_onshore: available_area
    #      [ ] power_densities
    # [ ] awe_shallow_fw1, wind_offshore: energy_cap_max group constraint
    #      [x] power densities: 8MW/km2
    # [ ] awe_deep_fw1, wind_floating: energy_cap_max group constraint
    #      [x] power densities: 8MW/km2

    # wind_onshore_monopoly, rooftop_pv: energy_cap_max
    df_energy_cap_max = get_energy_cap(model)
    print("df_energy_cap_max", df_energy_cap_max)
    # power_densities

    # open_field_pv, wind_onshore_competing, awe_onshore: available_area
    df_available_area_onshore = get_available_area_onshore(model)
    print("df_available_area_onshore", df_available_area_onshore)

    # power_densities

    # awe_shallow_fw1, wind_offshore: energy_cap_max group constraint
    df_offshore = get_offshore(model)
    print("df_offshore", df_offshore.index)
    joined = df_energy_cap_max.join(df_available_area_onshore).join(df_offshore)

    joined.to_csv(snakemake.output[0])
