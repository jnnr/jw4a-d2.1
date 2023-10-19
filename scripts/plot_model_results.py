import pandas as pd
import plotnine as pn


if __name__ == "__main__":
    energy_cap = pd.read_csv(snakemake.input[0])
    energy_cap["locs"] = energy_cap["locs"].replace({"_\\d+$": ""}, regex=True)
    energy_cap.groupby(["locs", "techs"]).sum()

    renewables = ["open_field_pv", "roof_mounted_pv", "wind", "awe"]
    (
        pn.ggplot(energy_cap.loc[energy_cap.techs.apply(lambda x: any([ren in x for ren in renewables]), 1),:])
        + pn.geom_col(pn.aes(x="locs", y="energy_cap", fill="techs"))
        + pn.labs(x="Region", y="Capacity")
        + pn.scale_color_discrete(guide=False)
    ).save(snakemake.output[0], dpi=300, height=6, width=12, facecolor="w", transparent=False)
