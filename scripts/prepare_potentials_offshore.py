r"""
Inputs
------
- offshore boundaries
- natura2000 exclusion
- other exclusion criteria

Outputs
-------
- potentials-wind-offshore-deep.csv
- potentials-wind-offshore-shallow.csv
- capacityfactors-wind-offshore-deep.csv
- capacityfactors-wind-offshore-shallow.csv
- shapefiles for map (nice to have, not straight forward with atlite, map can be done easily though

Description
-----------
Prepare offshore potentials in these steps:
"""
import geopandas as gpd
from pathlib import Path
import pandas as pd

if __name__ == "__main__":
    if "snakemake" not in globals():
        from lib.helpers import mock_snakemake

        snakemake = mock_snakemake("prepare_potentials_offshore")

    MW_TO_100GW = 1e-5 # (100GW/MW)
    M2_TO_KM2 = 1e-6 # (km2/m2)

    areas = pd.read_csv(snakemake.input.areas, index_col=0)

    potentials = areas.loc[:, ["available_area"]]

    potentials["energy_cap_max"] = potentials["available_area"] * snakemake.config["power_density"]
    potentials["energy_cap_max"] = potentials["energy_cap_max"] * M2_TO_KM2 * MW_TO_100GW

    # map eez to eurospores regions
    mapping = pd.read_csv(snakemake.input.mapping, index_col=0)
    mapping = mapping["country_code"].to_dict()
    potentials = potentials.loc[mapping.keys(), :]
    potentials = potentials.rename(index=mapping)

    potentials.to_csv(snakemake.output.potentials)
