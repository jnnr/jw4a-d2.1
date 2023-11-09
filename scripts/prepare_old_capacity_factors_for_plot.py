import pandas as pd

if __name__ == "__main__":
    if "snakemake" not in globals():
        from lib.helpers import mock_snakemake

        snakemake = mock_snakemake("prepare_old_capacity_factors_for_plot")

    capfac_existing = pd.read_csv(snakemake.input[0], index_col=0, parse_dates=True)

    capfac_existing.columns.name = "id"
    capfac_existing.index.name = "time"
    capfac_existing = capfac_existing.stack("id")
    capfac_existing.name = "__xarray_dataarray_variable__"
    capfac_existing = capfac_existing.to_xarray()
    capfac_existing = capfac_existing.to_dataset()

    capfac_existing.to_netcdf(snakemake.output[0])
