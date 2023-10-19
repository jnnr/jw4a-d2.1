r"""
Stripped down from pypsa-eur: https://github.com/PyPSA/pypsa-eur/blob/master/scripts/build_cutout.py

"""
import atlite
import pandas as pd


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_cutout")

    cutout_params = snakemake.config["build_cutout"]["cutout_params"]

    snapshots = pd.date_range(freq="h", **snakemake.config["snapshots"])
    time = [snapshots[0], snapshots[-1]]
    cutout_params["time"] = slice(*cutout_params.get("time", time))

    cutout_params["x"] = slice(*cutout_params["x"])
    cutout_params["y"] = slice(*cutout_params["y"])

    features = cutout_params.pop("features", None)
    cutout = atlite.Cutout(snakemake.output[0], **cutout_params)
    cutout.prepare(features=features)
