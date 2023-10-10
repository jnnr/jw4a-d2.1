import geopandas as gpd
import pandas as pd


def get_mapping_advanced(eez, eurospores):
    raise UserWarning("Not tested and probably not working.")
    # more sophisticated mapping of eez to eurospores. Problem: Spanish EEZ's geometry includes both atlantic and mediterranean sea, should be split (but before capacityfactors are calculated)
    # need to buffer the eurospores to make sure they intersect
    buffered_eurospores = eurospores.loc[:, ["id", "geometry"]]
    buffered_eurospores.geometry = buffered_eurospores.buffer(10000)
    intersecting = gpd.sjoin(eez.loc[eez["iso_sov2"].isna()], buffered_eurospores, how="left", predicate="intersects")

    # merge geometry back in
    mapping = intersecting.loc[:, ["id_left", "geometry", "id_right"]].merge(buffered_eurospores.loc[:, ["id", "geometry"]].rename(columns={"id": "id_right", "geometry": "geometry_right"}), on="id_right", how="left"
        )

    # rank by area of intersection
    mapping["rank"] = mapping.apply(lambda x: x["geometry"].intersection(x["geometry_right"]).area, axis=1)

    # keep those with highest overlap
    mapping = mapping.sort_values("rank", ascending=False).groupby("id_left").first()
    
    return mapping


def get_mapping_simple(eez, eurospores):
    # very simple mapping of eez to eurospores.
    mapping = eez.loc[eez["iso_sov2"].isna()][["id", "iso_sov1"]]
    mapping["iso_sov1"] = mapping["iso_sov1"].apply(lambda x: x + "_1")
    mapping = mapping.rename(columns={"iso_sov1": "country_code"})
    return mapping


if __name__ == "__main__":
    eez = gpd.read_file(snakemake.input.eez)
    eurospores = gpd.read_file(snakemake.input.eurospores)

    mapping = get_mapping_simple(eez, eurospores)

    mapping.to_csv(snakemake.output[0], index=False)
