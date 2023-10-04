rule draw_map:
    conda: "../envs/plot.yaml"
    input: 
        path_water_depth="data/potentials_offshore/gebco_2023_sub_ice_topo/GEBCO_2023_sub_ice_topo.nc",
        path_boundaries_eez="build/shapes/eez.geojson",
        path_boundaries_onshore="data/europe-98-zones.geojson/europe-98-zones.geojson"
    output: "build/plots/map.png"
    script: "../scripts/draw_map.py"

rule draw_plots_capacityfactors:
    conda: "../envs/plot.yaml"
    input:
        path_capacity_factors = "build/capacity_factors/capacity_factors_{tech}.nc",
        path_boundaries = "build/shapes/eez.geojson"
    output:
        path_plot = "build/plots/load_duration_wind_{tech}.png",
        path_plot_average = "build/plots/capacity_factor_average_{tech}.png"
    script: "../scripts/draw_plots_capacityfactors.py"
