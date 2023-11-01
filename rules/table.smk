rule table_area_potential_offshore:
    input:
        area_offshore_deep = "build/availability/area_offshore_deep.csv",
        area_offshore_shallow = "build/availability/area_offshore_shallow.csv"
    output: "build/tables/area_potential_offshore.csv"
    script: "../scripts/table_area_potential.py"

rule table_area_constraints:
    conda: "../envs/calliope.yaml"
    input: "run-prebuilt-sector-coupled-euro-calliope/build/eurospores/outputs/2016_res_6h_w_noveltech.nc"
    output: "build/tables/area_constraints.csv"
    script: "../scripts/table_area_constraints.py"
