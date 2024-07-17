#!/usr/bin/env python

import click
import geopandas as gpd
import numpy as np
import rasterio as rio

from illum import MultiScaleData as MSD


@click.command(name="convert")
@click.option("--vector/--raster", "-v/-r", default=True, help="Output type.")
@click.option("-log", is_flag=True, help="Logarithmic scale (base 10)")
@click.option("-area", is_flag=True, help="Normalized by pixel area (in km²)")
@click.argument("filename", type=click.Path(exists=True))
@click.argument("outname")
def CLI_convert(filename, outname, vector, log, area):
    """Convert an Illumina HDF file to a georeferenced format.

    Converts FILENAME to OUTNAME.EXT where ext is defined based on the output
    format. The output format is either vector (GeoJSON) or raster (Tiff).
    """
    convert(filename, outname, vector, log, area)


def convert(filename, outname, vector=True, log=False, area=False):
    hdf = MSD.Open(filename)
    hdf.set_buffer(-1)
    hdf.set_overlap(-1)

    if not vector:
        profile = rio.default_gtiff_profile.copy()

        for i, data in enumerate(hdf):
            pix_size = hdf._attrs["layers"][i]["pixel_size"]
            b = hdf._attrs["layers"][i]["buffer"]
            xmin = hdf._attrs["layers"][i]["xmin"] + b * pix_size
            ymax = hdf._attrs["layers"][i]["ymax"] - b * pix_size
            if b != 0:
                data = data[b:-b, b:-b]
            if area:
                data /= (hdf.pixel_size(i) / 1000.0) ** 2
            if log:
                data = np.log10(data)

            profile = dict(
                driver="GTiff",
                dtype=rio.float64,
                nodata=-np.inf if log else -1,
                width=data.shape[1],
                height=data.shape[0],
                count=1,
                crs=rio.CRS.from_user_input(hdf._attrs['srs']),
                transform=rio.transform.from_origin(xmin,ymax,pix_size,pix_size),
            )

            with rio.open(f"{outname}_{i}.tif", 'w', **profile) as dst:
                dst.write(data.astype(rio.float64), 1)

    elif vector:
        points = {"x": [], "y": [], "val": []}
        for i, data in enumerate(hdf):
            xmin = hdf._attrs["layers"][i]["xmin"]
            ymax = hdf._attrs["layers"][i]["ymax"]
            pix_size = hdf._attrs["layers"][i]["pixel_size"]

            pts = np.where(hdf[i] != -1)
            points["x"].extend((pts[1] + 0.5) * pix_size + xmin)
            points["y"].extend(ymax - (pts[0] + 0.5) * pix_size)
            data = hdf[i][hdf[i] != -1]
            if area:
                data /= (hdf.pixel_size(i) / 1000.0) ** 2
            if log:
                data = np.log10(data)
            points["val"].extend(data)

        gdf = gpd.GeoDataFrame(
            points,
            crs=hdf._attrs["srs"],
            geometry=gpd.points_from_xy(points["x"], points["y"]),
        )

        gdf.to_file(outname + ".geojson", driver="GeoJSON")
