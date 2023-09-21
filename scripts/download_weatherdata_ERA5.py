r"""
Downloads wind timeseries from ERA5
"""
import cdsapi
import pandas as pd
import os


def download_ERA5_ml(lon_west, lon_east, lat_north, lat_south, filepath, client):
    print(f"Downloading from ERA5, model level. Saving to {filepath}")
    print(f"lon_west={lon_west}, lon_east={lon_east}, lat_north={lat_north}, lat_south={lat_south}")
    client.retrieve(
        "reanalysis-era5-complete",
        {
            "date": "2013-01-01/to/2018-12-31",
            "levelist": "127",
            "levtype": "ml",
            "param": "131/132",
            "stream": "oper",
            "time": "00:00:00/01:00:00/02:00:00/03:00:00/04:00:00/05:00:00/06:00:00/07:00:00/08:00:00/09:00:00/10:00:00/11:00:00/12:00:00/13:00:00/14:00:00/15:00:00/16:00:00/17:00:00/18:00:00/19:00:00/20:00:00/21:00:00/22:00:00/23:00:00",
            "type": "an",
            "area": [lat_north, lon_west, lat_south, lon_east],
            "grid": "0.25/0.25",
            "format": "netcdf",
        },
        filepath,
    )


def download_ERA5_pl(lon_west, lon_east, lat_north, lat_south, filepath, client):
    print("Downloading from ERA5, pressure level.")
    print(lon_west, lon_east, lat_north, lat_south, filepath, client)
    client.retrieve(
        "reanalysis-era5-pressure-levels",
        {
            "product_type": "reanalysis",
            "variable": [
                "u_component_of_wind",
                "v_component_of_wind",
            ],
            "pressure_level": "975",
            "year": [
                "2013",
                "2014",
                "2015",
                "2016",
                "2017",
                "2018",
            ],
            "month": [
                "01",
                "02",
                "03",
                "04",
                "05",
                "06",
                "07",
                "08",
                "09",
                "10",
                "11",
                "12",
            ],
            "day": [
                "01",
                "02",
                "03",
                "04",
                "05",
                "06",
                "07",
                "08",
                "09",
                "10",
                "11",
                "12",
                "13",
                "14",
                "15",
                "16",
                "17",
                "18",
                "19",
                "20",
                "21",
                "22",
                "23",
                "24",
                "25",
                "26",
                "27",
                "28",
                "29",
                "30",
                "31",
            ],
            "time": [
                "00:00",
                "01:00",
                "02:00",
                "03:00",
                "04:00",
                "05:00",
                "06:00",
                "07:00",
                "08:00",
                "09:00",
                "10:00",
                "11:00",
                "12:00",
                "13:00",
                "14:00",
                "15:00",
                "16:00",
                "17:00",
                "18:00",
                "19:00",
                "20:00",
                "21:00",
                "22:00",
                "23:00",
            ],
            "area": [lat_north, lon_west, lat_south, lon_east],
            "format": "netcdf",
        },
        filepath,
    )


if __name__ == "__main__":
    target_dir = snakemake.output.target_dir

    cdsclient = cdsapi.Client(timeout=600, quiet=False, debug=True)

    download_ERA5_ml(-17., 37., 72., 32., target_dir, cdsclient)