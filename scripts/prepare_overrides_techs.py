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
    
Example :file:`advanced-wind-techs.yaml`:

.. code-block:: yaml
    techs:
        awe_shallow_fw1: 
        essentials:
            name: AWE shallow fw1
            color:
            parent: supply
            carrier: electricity
        constraints: 
            resource: file=AWE_shallow_fw1.csv
            resource_unit: energy_per_cap
            lifetime: 30
        costs:
            monetary:
                energy_cap: 153.35 #ampyx
                om_annual: 4.79
                om_prod: 0.00024

Example :file:`advanced-wind-locations.yaml` for setting `energy_cap_max`:

.. code-block:: yaml
    locations:
        ALB_1: # Albania
            techs:
                wind_onshore_monopoly:
                    constraints:
                        energy_cap_max: 0.6879766264903546 # (100,000 MW)

                        
Example :file:`advanced-wind-locations.yaml` for setting group constraints:  

.. code-block:: yaml                 
    overrides:
        offshorewind_cap_max:
            group_constraints:
                offshorewind_BEL:
                    techs: [wind_offshore, awe_shallow_fw1, awe_shallow_fw2]
                    locs: [BEL]
                    energy_cap_max: 0.11  # (100,000 MW)
"""
from pathlib import Path
from template import parametrise_template


if __name__ == "__main__":
    if "snakemake" not in globals():
        template_path = Path("templates")
        output_path = Path("build")
        templates = [
            "advanced-wind-techs.yaml"
        ]

    # open scalar data

    # parametrize template
    for template in templates:
        parametrise_template(template_path / template, output_path / template)
