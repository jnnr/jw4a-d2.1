rule draw_map:
    input: 
        path_water_depth="data/potentials_offshore/gebco_2023_sub_ice_topo/GEBCO_2023_sub_ice_topo.nc",
        path_boundaries_eez="build/shapes/eez.geojson",
        path_boundaries_onshore="data/europe-98-zones.geojson/europe-98-zones.geojson"
    output: "build/plots/map.png"
    script: "../scripts/draw_map.py"