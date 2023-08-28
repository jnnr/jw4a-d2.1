rule prepare_coordinates_from_shapefile:
    input: 
        path_geoboundaries="data/europe-98-zones.geojson/europe-98-zones.geojson"
    output: 
        path_coordinates="build/coordinates/europe-98-zones.csv"
    script: "../scripts/prepare_coordinates.py"

rule prepare_offshore_boundaries:
    input:
        eez="data/EEZ_land_union_v3_202003/EEZ_land_union_v3_202003",
    output:
        boundaries_offshore="build/shapes/boundaries_offshore.geojson"
    script: "../scripts/prepare_offshore_boundaries.py"