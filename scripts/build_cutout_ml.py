import xarray as xr
import numpy as np


def _rename_and_clean_coords(ds, add_lon_lat=True):
    """
    Rename 'longitude' and 'latitude' columns to 'x' and 'y' and fix roundings.

    Optionally (add_lon_lat, default:True) preserves latitude and
    longitude columns as 'lat' and 'lon'.
    """
    ds = ds.rename({"longitude": "x", "latitude": "y"})
    # round coords since cds coords are float32 which would lead to mismatches
    ds = ds.assign_coords(
        x=np.round(ds.x.astype(float), 5), y=np.round(ds.y.astype(float), 5)
    )

    if add_lon_lat:
        ds = ds.assign_coords(lon=ds.coords["x"], lat=ds.coords["y"])

    return ds


if __name__ == "__main__":
    modellevel = xr.open_dataset(snakemake.input)

    modellevel.attrs["module"] = "era5"
    modellevel.attrs["dx"] = 0.25
    modellevel.attrs["dy"] = 0.25
    modellevel.attrs["prepared_features"] = ["u", "v"]
    modellevel = _rename_and_clean_coords(modellevel)

    for v in modellevel:
        modellevel[v].attrs["feature"] = v
        modellevel[v].attrs["module"] = "era5-ml"

    modellevel["wind_speed"] = np.sqrt(modellevel["u"] ** 2 + modellevel["v"] ** 2).assign_attrs(
        units=modellevel["u"].attrs["units"], long_name="Wind speed"
    )

    modellevel.to_netcdf(snakemake.output)
    modellevel.close()
