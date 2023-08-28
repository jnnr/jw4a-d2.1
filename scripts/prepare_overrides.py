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