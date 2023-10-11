rule draw_map:
    conda: "../envs/plot.yaml"
    input: 
        water_depth="data/potentials_offshore/gebco_2023_sub_ice_topo/GEBCO_2023_sub_ice_topo.nc",
        boundaries_eez="build/shapes/eez.geojson",
        boundaries_onshore="data/europe-98-zones.geojson/europe-98-zones.geojson",
        natura2000="data/potentials_offshore/natura2000_areas/eea_v_3035_100_k_natura2000_p_2021_v12_r01/SHP/Natura2000_end2021_rev1_epsg3035.shp",
        cutout_era5="build/cutouts/cutout-era5.nc",
        cutout_era5_model_level="build/cutouts/cutout-era5-model-level_adapted.nc"
    output: 
        areas="build/plots/map.png",
        wind_speeds_era5="build/plots/map_wind_speeds_era5.png",
        wind_speeds_era5_model_level="build/plots/map_wind_speeds_model_level.png"
    script: "../scripts/draw_map.py"

def get_path_boundaries(wildcards):
    path_boundaries = config["draw_plots_capacityfactors"][wildcards.tech]
    return path_boundaries

rule draw_plots_capacityfactors:
    conda: "../envs/plot.yaml"
    input:
        path_capacity_factors = "build/capacity_factors/capacity_factors_{tech}.nc",
        path_boundaries = get_path_boundaries
    output:
        path_plot = "build/plots/load_duration_wind_{tech}.png",
        path_plot_average = "build/plots/capacity_factor_average_{tech}.png"
    script: "../scripts/draw_plots_capacityfactors.py"


rule prepare_old_capacity_factors_for_plot:
    input: "run-prebuilt-sector-coupled-euro-calliope/build/pre-built/model/eurospores/capacityfactors-{tech}.csv"
    output: "build/capacity_factors/capacity_factors_old_{tech}.nc"
    wildcard_constraints:
        tech = "wind-offshore|wind-onshore"
    script: "../scripts/prepare_old_capacity_factors_for_plot.py"

