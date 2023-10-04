import atlite
import functools
import numpy as np
import geopandas as gpd
from pathlib import Path
import pandas as pd


def compute_shape_availabilities(geometries, excluder):
    assert isinstance(geometries, gpd.GeoDataFrame)
    availabilities = pd.DataFrame(geometries.drop(columns="geometry"))
    for id in geometries.index:
        geometry = geometries.loc[[id]]

        masked, _ = excluder.compute_shape_availability(geometry)
        availability = masked.sum() * excluder.res**2

        area = geometry.geometry.to_crs(excluder.crs).item().area

        available_share = availability / area
        availabilities.loc[id, ["area", "available_area", "available_share"]] = [area, availability, available_share]

    return availabilities


if __name__ == "__main__":
    if "snakemake" not in globals():
        from lib.helpers import mock_snakemake

        snakemake = mock_snakemake("build_availabilitymatrix")

    # prepare boundaries
    boundaries = gpd.read_file(snakemake.input.boundaries)

    # prepare cutout
    cutout = atlite.Cutout(snakemake.input.cutout)

    # prepare excluder
    excluder_shallow = atlite.ExclusionContainer()
    excluder_shallow.add_geometry(snakemake.input.natura2000)
    max_depth = 60
    func = functools.partial(np.greater, -max_depth)
    excluder_shallow.add_raster(snakemake.input.gebco, crs=4326, codes=func)

    excluder_deep = atlite.ExclusionContainer()
    excluder_deep.add_geometry(snakemake.input.natura2000)
    min_depth = 60
    func = functools.partial(np.less, -max_depth)
    excluder_deep.add_raster(snakemake.input.gebco, crs=4326, codes=func)

    # prepare availability matrix
    if not Path(snakemake.output.area_deep).parent.exists():
        Path(snakemake.output.area_deep).parent.mkdir(parents=True)

    availability = cutout.availabilitymatrix(boundaries, excluder_deep)
    availability.to_netcdf(snakemake.output.availability_deep)
    area_deep = compute_shape_availabilities(boundaries, excluder_deep)
    area_deep.to_csv(snakemake.output.area_deep)
    

    availability = cutout.availabilitymatrix(boundaries, excluder_shallow)
    availability.to_netcdf(snakemake.output.availability_shallow)
    area_shallow = compute_shape_availabilities(boundaries, excluder_shallow)
    area_shallow.to_csv(snakemake.output.area_shallow)
