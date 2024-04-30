rule build_cutout_model_level:
    input: "build/cutouts/cutout-era5-model-level.nc"
    output: "build/cutouts/cutout-era5-model-level_adapted.nc"
    script: "../scripts/build_cutout_model_level.py"

rule build_availabilitymatrix:
    input:
        boundaries = "build/shapes/eez.geojson",
        cutout = "build/cutouts/cutout-era5-model-level.nc",
        natura2000 = "data/potentials_offshore/natura2000_areas/eea_v_3035_100_k_natura2000_p_2021_v12_r01/SHP/Natura2000_end2021_rev1_epsg3035.shp",
        shipdensity = "data/potentials_offshore/shipdensity_global/shipdensity_global.tif",
        gebco = "data/potentials_offshore/gebco_2023_sub_ice_topo/GEBCO_2023_sub_ice_topo.nc"
    output:
        availability_deep = "build/availability/availability_offshore_deep.nc",
        availability_shallow = "build/availability/availability_offshore_shallow.nc",
        area_deep = "build/availability/area_offshore_deep.csv",
        area_shallow = "build/availability/area_offshore_shallow.csv"
    script: "../scripts/build_availabilitymatrix.py"

def get_inputs(wildcards):
    inputs = config["prepare_capacity_factors"][wildcards.tech]
    return inputs

rule prepare_capacity_factors:
    input: 
        unpack(get_inputs)
    output: "build/capacity_factors/capacity_factors_{tech}.nc"
    wildcard_constraints:
        tech = "|".join(config["prepare_capacity_factors"].keys())
    script: "../scripts/prepare_capacity_factors.py"
