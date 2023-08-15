rule download_ERA5_weather_data:
    input: 
        path_coordinates="data/coordinates/coordinates.csv"
    output: 
        target_dir=directory("build/weather_data")
    script: "../scripts/download_weatherdata_ERA5.py"

rule prepare_capacity_factors:
    input: 
        windspeeds="build/weather_data",
        powercurve="data/powercurve.csv"
    output: 
        capacity_factors="build/capacity_factors"
    script: "../scripts/prepare_capacity_factors.py"
