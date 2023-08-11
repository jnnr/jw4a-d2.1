rule prepare_coordinates_from_shapefile:
    input: 
        path_geoboundaries="data/europe-98-zones.geojson/europe-98-zones.geojson"
    output: 
        path_coordinates="build/coordinates/europe-98-zones"
    script: "../scripts/prepare_coordinates.py"