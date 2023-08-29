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

    # boundaries_offshore = gpd.read_file(snakemake.input.boundaries_offshore)
    # water_depth = gpd.read_file(snakemake.input.water_depth)
    # natura2000 = gpd.read_file(snakemake.input.natura2000)

    # TODO: Why not one df?
    potentials_offshore_deep = pd.DataFrame()
    potentials_offshore_shallow = pd.DataFrame()

    potentials_offshore_deep.to_csv(snakemake.output.potentials_offshore_deep)
    potentials_offshore_shallow.to_csv(snakemake.output.potentials_offshore_shallow)