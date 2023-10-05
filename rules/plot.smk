rule draw_map:
    conda: "../envs/plot.yaml"
    input: 
        path_water_depth="data/potentials_offshore/gebco_2023_sub_ice_topo/GEBCO_2023_sub_ice_topo.nc",
        path_boundaries_eez="build/shapes/eez.geojson",
        path_boundaries_onshore="data/europe-98-zones.geojson/europe-98-zones.geojson",
        path_natura2000 = "data/potentials_offshore/natura2000_areas/eea_v_3035_100_k_natura2000_p_2021_v12_r01/SHP/Natura2000_end2021_rev1_epsg3035.shp"
    output: "build/plots/map.png"
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
