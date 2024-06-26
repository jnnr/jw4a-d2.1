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

rule draw_boxplot_capacityfactors:
    conda: "../envs/plot.yaml"
    input:
        offshore_deep_awe="build/capacity_factors/capacity_factors_offshore_deep_awe.nc",
        offshore_shallow_awe="build/capacity_factors/capacity_factors_offshore_shallow_awe.nc",
        old_wind_offshore="build/capacity_factors/capacity_factors_old_wind-offshore.nc",
        old_wind_onshore="build/capacity_factors/capacity_factors_old_wind-onshore.nc",
        onshore_awe="build/capacity_factors/capacity_factors_onshore_awe.nc"
    output:
        path_plot = "build/plots/boxplot_capacityfactors.png",
    script: "../scripts/draw_boxplots_capacityfactors.py"

rule prepare_old_capacity_factors_for_plot:
    input: "run-prebuilt-sector-coupled-euro-calliope/build/pre-built/model/eurospores/capacityfactors-{tech}.csv"
    output: "build/capacity_factors/capacity_factors_old_{tech}.nc"
    wildcard_constraints:
        tech = "wind-offshore|wind-onshore"
    script: "../scripts/prepare_old_capacity_factors_for_plot.py"

rule plot_capacity_factor_distribution:
    input:
        offshore_deep_awe="build/capacity_factors/capacity_factors_offshore_deep_awe.nc",
        offshore_shallow_awe="build/capacity_factors/capacity_factors_offshore_shallow_awe.nc",
        old_wind_offshore="build/capacity_factors/capacity_factors_old_wind-offshore.nc",
        old_wind_onshore="build/capacity_factors/capacity_factors_old_wind-onshore.nc",
        onshore_awe="build/capacity_factors/capacity_factors_onshore_awe.nc",
    output:
        histogram="build/plots/capacity_factors_histogram.png",
        regional_histogram="build/plots/capacity_factors_regional_histogram.png",
        # boxplot="build/plots/capacity_factors_boxplot.png",
        regional_boxplot="build/plots/capacity_factors_regional_boxplot.png",
        description="build/plots/capacity_factors_description.csv"
    script: "../scripts/plot_capacity_factor_distribution.py"

rule postprocess_model_results:
    conda: "../envs/calliope.yaml"
    input: ancient("run-prebuilt-sector-coupled-euro-calliope/build/eurospores/outputs/{year}_{resolution}_{model_resolution}_{scenario}.nc")
    output:
        energy_cap="build/postprocessed_results/{year}_{resolution}_{model_resolution}_{scenario}_energy_cap.csv",
        tech_names="build/postprocessed_results/{year}_{resolution}_{model_resolution}_{scenario}_tech_names.json"
    script: "../scripts/postprocess_model_results.py"


rule plot_model_results:
    conda: "../envs/plot.yaml"
    input:
        energy_cap="build/postprocessed_results/{year}_{resolution}_{model_resolution}_{scenario}_energy_cap.csv",
        tech_names="build/postprocessed_results/{year}_{resolution}_{model_resolution}_{scenario}_tech_names.json"
    output: "build/plots/{year}_{resolution}_{model_resolution}_{scenario}_energy_cap.png"
    script: "../scripts/plot_model_results.py"