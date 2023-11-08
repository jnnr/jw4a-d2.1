rule table_area_potential_offshore:
    input:
        area_offshore_deep = "build/availability/area_offshore_deep.csv",
        area_offshore_shallow = "build/availability/area_offshore_shallow.csv"
    output: "build/tables/area_potential_offshore.csv"
    script: "../scripts/table_area_potential.py"

rule table_area_constraints:
    conda: "../envs/calliope.yaml"
    input: ancient("run-prebuilt-sector-coupled-euro-calliope/build/eurospores/outputs/{year}_{resolution}_{model_resolution}_{scenario}.nc")
    output: "build/tables/{year}_{resolution}_{model_resolution}_{scenario}_area_constraints.csv"
    script: "../scripts/table_area_constraints.py"
