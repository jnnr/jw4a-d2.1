import pandas as pd
import numpy as np


def format_tuple(tupl):
    relevant_entries = [item for item in list(tupl) if isinstance(item, str)]
    result = "/".join(relevant_entries)
    return result


if __name__ == "__main__":

    def prepare_table(path_area_potential):
        METER_PER_KM = 1e6
        area_potential = pd.read_csv(path_area_potential)

        table = area_potential.loc[:, ["area_km2", "available_area", "available_share"]]

        table["shortname"] = area_potential[["iso_sov1", "iso_sov2"]].apply(format_tuple, axis=1)

        # format numbers
        table.loc[:, "available_share"] = np.round(table.loc[:, "available_share"], 2)
        table.loc[:, "available_area"] = table.loc[:, "available_area"].apply(lambda x: f"{x/METER_PER_KM:.0f}")

        table = table.set_index(["shortname", "area_km2"])

        return table
    
    table_deep = prepare_table(snakemake.input.area_offshore_deep)
    table_shallow = prepare_table(snakemake.input.area_offshore_shallow)

    table_joined = table_deep.join(table_shallow, lsuffix="_deep", rsuffix="_shallow")

    table_joined = table_joined.sort_values(by=["area_km2"], ascending=False)

    table_joined.columns = [col.replace("_", " ") for col in table_joined.columns]

    table_joined.to_csv(str(snakemake.output))
