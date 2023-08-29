rule prepare_potentials_offshore:
    input:
        offshore_boundaries="build/shapes/boundaries_offshore.geojson",
        natura2000="data/potentials_offshore/natura2000_areas",
        water_depth="data/potentials_offshore/gebco_2023_sub_ice_topo"
    output:
        potentials_offshore_deep="build/potentials_offshore/potentials_offshore_deep.csv",
        potentials_offshore_shallow="build/potentials_offshore/potentials_offshore_shallow.csv"
    script: "../scripts/prepare_potentials_offshore.py"