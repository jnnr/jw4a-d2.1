# coding: utf-8
r"""
Inputs
------
- cost energy_cap
- cost interest_rate
- lifetime
- regions
- Capacity factor timeseries
- filepath to store yaml file

Outputs
-------
- advanced-wind-techs.yaml file containing tech.
- advanced-wind-locations.yaml

Description
-------------
Create a yaml file that adds AWE to the model.

These techs are added to the model: 
Airborne wind onshore
Airborne wind shallow water
Airborne wind floating

We assume
- that offshore wind and awe shallow water share the same space using a group constraint (however, Hidde assumes shared energy_cap_max, not area!)
- that offshore floating and awe deep water share the same space (however, Hidde assumes shared energy_cap_max, not area!)
"""
import sys
import pandas as pd
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent))
from lib.template import parametrise_template
import xarray as xr


def prepare_capacity_factors(capfac, destination, name, rename_columns=None):
    if not destination.exists():
        destination.mkdir(parents=True)

    # chunk the data into years and save it as csv
    years = list(set(capfac.time.dt.year.values))
    for year in years:
        capfac_year = capfac.sel(time=str(year))
        df_capfac_year = capfac_year.to_dataframe()["__xarray_dataarray_variable__"].unstack("id")

        if rename_columns is not None:
            df_capfac_year = df_capfac_year.rename(columns=rename_columns)

        df_capfac_year.to_csv(destination / f"{name}_{year}.csv")

def prepare_areas():
    pass


if __name__ == "__main__":
    if "snakemake" not in globals():
        from lib.helpers import mock_snakemake

        snakemake = mock_snakemake("prepare_overrides")

    # Generate output path
    path_output = Path(snakemake.output[0])

    if not path_output.exists():
        path_output.mkdir(parents=True)

    # parametrize template for locations
    for template, path_potentials in zip(snakemake.input.templates_locations, snakemake.input.potentials):
        path_target = path_output / Path(template).name
        potentials = pd.read_csv(path_potentials, index_col=0)
        parametrise_template(template, path_target, potentials=potentials)

    # parametrize template for techs
    parametrise_template(snakemake.input.template_techs, path_output / Path(snakemake.input.template_techs).name)


    for path_capacity_factors in snakemake.input.capacity_factors:
        path_capacity_factors = Path(path_capacity_factors)
        capacity_factors = xr.load_dataset(path_capacity_factors)

        mapping = None
        if "offshore" in str(path_capacity_factors):
            mapping = pd.read_csv(snakemake.input.mapping, index_col=0)
            mapping.index = mapping.index.astype(str)
            mapping = mapping["country_code"].to_dict()
    
        prepare_capacity_factors(capacity_factors, destination=path_output / "timeseries", name=path_capacity_factors.stem, rename_columns=mapping)
    