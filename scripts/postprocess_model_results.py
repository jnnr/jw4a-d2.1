import json

import calliope
import pandas as pd


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


def get_tech_names(model):
    return model.inputs.names.to_pandas().to_dict()


if __name__ == "__main__":
    if "snakemake" not in globals():
        from lib.helpers import mock_snakemake

        snakemake = mock_snakemake("postprocess_model_results")

    model = calliope.read_netcdf(snakemake.input[0])

    energy_cap = get_tidy_data(model, "energy_cap").sort_values(by=["locs", "techs"])

    energy_cap_max = get_tidy_data(model, "energy_cap_max").sort_values(
        by=["locs", "techs"]
    )

    df = energy_cap.merge(energy_cap_max, on=["locs", "techs"])

    df.to_csv(snakemake.output.energy_cap, index=False)

    tech_names = get_tech_names(model)
    with open(snakemake.output.tech_names, "w") as f:
        json.dump(tech_names, f)
