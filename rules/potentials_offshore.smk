rule prepare_potentials_offshore:
    input:
        areas="build/availability/area_offshore_{depth}.csv",
        mapping="build/shapes/map_eez_eurospores.csv"
    output:
        potentials="build/potentials_offshore/potentials_offshore_{depth}.csv"
    wildcard_constraints:
        depth="deep|shallow"
    script: "../scripts/prepare_potentials_offshore.py"
