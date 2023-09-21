rule build_cutout:
    input: "build/cutouts/cutout-era5-model-level.nc"
    output: "build/cutouts/cutout-era5-model-level_adapted.nc"
    script: "../scripts/build_cutout_ml.py"

rule prepare_capacity_factors:
    input: 
        windspeed="build/weather_data",
        powercurve="data/power_curves/AWE_500kw_softwing.csv"
    output: 
        capacity_factors="build/capacity_factors/AWE_500kw_softwing.csv"
    script: "../scripts/prepare_capacity_factors.py"
