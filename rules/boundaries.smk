rule prepare_coordinates_from_shapefile:
    input: 
        path_geoboundaries="data/europe-98-zones.geojson/europe-98-zones.geojson"
    output: 
        path_coordinates="build/coordinates/europe-98-zones.csv"
    script: "../scripts/prepare_coordinates.py"

rule eez:
    # Adapted from euro-calliope
    message: "Clip exclusive economic zones to study area."
    input: "data/shapes/eez.zip"
    output: "build/shapes/eez.geojson"
    params:
        bounds="{x_min},{y_min},{x_max},{y_max}".format(**config["scope"]["spatial"]["bounds"]),
        countries=",".join(["'{}'".format(country) for country in config["scope"]["spatial"]["countries"]]),
    # conda: "../envs/geo.yaml"
    shadow: "minimal"
    shell:
        """
        fio cat --bbox {params.bounds} "zip://{input}"\
        | fio filter "f.properties.territory1 in [{params.countries}]"\
        | fio collect > {output}
        """

rule prepare_offshore_boundaries:
    input:
        eez="data/potentials_offshore/eez"
    output:
        boundaries_offshore="build/shapes/boundaries_offshore.geojson"
    script: "../scripts/prepare_offshore_boundaries.py"