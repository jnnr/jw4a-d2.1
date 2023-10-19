rule prepare_overrides:
    input: 
        potentials=[
            "build/potentials_offshore/potentials_offshore_deep.csv",
            "build/potentials_offshore/potentials_offshore_shallow.csv",
            ],
        capacity_factors=[
            "build/capacity_factors/capacity_factors_offshore_deep_awe.nc",
            "build/capacity_factors/capacity_factors_offshore_shallow_awe.nc",
            "build/capacity_factors/capacity_factors_onshore_awe.nc",
        ],
        templates_locations=[
            "data/templates/locations-wind-offshore-deep.yaml",
            "data/templates/locations-wind-offshore-shallow.yaml",
        ],
        template_techs="data/templates/techs-novel-wind.yaml",
        mapping="build/shapes/map_eez_eurospores.csv"
    output:
        directory("build/overrides")
    script: "../scripts/prepare_overrides.py"
