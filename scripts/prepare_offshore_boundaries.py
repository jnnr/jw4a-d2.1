r"""
Inputs
------
- eez

Outputs
-------
- offshore boundaries

Description
-----------
Prepare boundaries for offshore regions.

Based on the eez.
Can be further segmented by some segmentation approach.
"""
from pathlib import Path

import geopandas as gpd

if __name__ == "__main__":
    if "snakemake" not in globals():
        from lib.helpers import mock_snakemake

        snakemake = mock_snakemake("prepare_shapes")

    eez = gpd.read_file(Path(snakemake.input.eez) / "EEZ_Land_v3_202030.shp")

    gdf = eez

    gdf.to_file(snakemake.output[0])
