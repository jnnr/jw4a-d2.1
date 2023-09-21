rule prepare_capacity_factors:
    input: 
        windspeed="build/weather_data",
        powercurve="data/power_curves/AWE_500kw_softwing.csv"
    output: 
        capacity_factors="build/capacity_factors/AWE_500kw_softwing.csv"
    script: "../scripts/prepare_capacity_factors.py"
