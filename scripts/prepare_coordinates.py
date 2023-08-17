import geopandas as gpd
import pandas as pd
import numpy as np


def round_to(x, base=0.25):
    return np.round(x/base)*base


def get_centroids_from_gdf(gdf):
    # TODO fix this
    # UserWarning: Geometry is in a geographic CRS. Results from 'centroid' are likely incorrect.
    # Use 'GeoSeries.to_crs()' to re-project geometries to a projected CRS before this operation.

    # Produce a gdf in this format:
    # country, location, lon_west, lon_east, lat_north, lat_south
    coordinates = gpd.GeoDataFrame()
    coordinates["country"] = gdf["name"]
    coordinates["location"] = gdf["id"]

    centroids = gdf["geometry"].centroid

    coordinates["lon_west"] = round_to(centroids.x)
    coordinates["lat_north"] = round_to(centroids.y)
    coordinates["lon_east"] = coordinates["lon_west"] + 0.125
    coordinates["lat_south"] = coordinates["lat_north"] - 0.125

    return coordinates


if __name__ == "__main__":
    path_geoboundaries = snakemake.input.path_geoboundaries
    path_coordinates = snakemake.output.path_coordinates
    file_geoboundaries = open(path_geoboundaries)
    geoboundaries = gpd.read_file(file_geoboundaries)
    centroids = get_centroids_from_gdf(geoboundaries)
    centroids.to_csv(path_coordinates, index=False)
