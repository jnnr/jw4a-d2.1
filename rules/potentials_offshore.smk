rule prepare_potentials_offshore:
    input:
        areas="build/availability/area_offshore_{depth}.csv",
        mapping="build/shapes/map_eez_eurospores.csv"
    output:
        potentials="build/potentials_offshore/potentials_offshore_{depth}.csv"
    wildcard_constraints:
        depth="deep|shallow"
    script: "../scripts/prepare_potentials_offshore.py"


rule prepare_locations_onshore:
    input: "data/europe-98-zones.geojson/europe-98-zones.geojson"
    output: "build/locations_onshore/locations_onshore.csv"
    run:
        import geopandas as gpd
        regions_onshore = gpd.read_file(input[0]).loc[:, "id"]
        regions_onshore.to_csv(output[0], index=False)