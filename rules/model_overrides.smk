rule prepare_overrides:
    input: 
        potentials=[
            "data/potentials_dummy/potentials-wind-offshore-deep.csv",
            "data/potentials_dummy/potentials-wind-offshore-shallow.csv",
            ],
        templates_locations=[
            "data/templates/locations-wind-offshore-deep.yaml",
            "data/templates/locations-wind-offshore-shallow.yaml",
        ],
        template_techs="data/templates/techs-novel-wind.yaml",
    output:
        directory("build/overrides")
    script: "../scripts/prepare_overrides.py"
