rule download_ERA5_weather_data:
    input: 
        path_coordinates="build/coordinates/europe-98-zones.csv"
    output: 
        target_dir=directory("build/weather_data")
    script: "../scripts/download_weatherdata_ERA5.py"

rule prepare_capacity_factors:
    input: 
        windspeed="build/weather_data",
        powercurve="data/power_curves/AWE_500kw_softwing.csv"
    output: 
        capacity_factors="build/capacity_factors/AWE_500kw_softwing.csv"
    script: "../scripts/prepare_capacity_factors.py"
