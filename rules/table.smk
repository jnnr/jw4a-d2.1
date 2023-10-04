rule table_area_potential_offshore:
    input:
        area_offshore_deep = "build/availability/area_offshore_deep.csv",
        area_offshore_shallow = "build/availability/area_offshore_shallow.csv"
    output: "build/tables/area_potential_offshore.csv"
    script: "../scripts/table_area_potential.py"
