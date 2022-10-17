#!/usr/bin/env python3

import os
from math import ceil

import numpy as np
import pyproj
import rasterio.features
import shapely.geometry
import shapely.ops
from progressbar import progressbar

import illum.utils as u


def roads(domain, res=1, xsize=10000, ysize=10000, buffer=200, epsg=None):
    print("Analyzing domain.")
    rst = rasterio.open(domain)
    coords = rasterio.transform.xy(rst.transform, *np.where(rst.read(1)))
    if epsg is None:
        epsg = u.estimate_utm_epsg(*coords)
    print(f"Using epsg:{epsg} projection.")
    bounds = shapely.geometry.MultiPoint(
        tuple(zip(*u.transform(t_crs=epsg)(*coords)))
    ).minimum_rotated_rectangle

    print("Splitting domain.")
    res = 1
    buffer = ceil(buffer / res) * res
    tiles, shape = u.fishnet(bounds, res=res, xsize=xsize, ysize=ysize)
    ny, nx = shape
    nb = int(buffer // res)
    mask = slice(nb, -nb), slice(nb, -nb)
    Ty, Tx = np.max(list(tiles.keys()), 0) + 1
    print(f"Domain split in {len(tiles)} ({Ty}x{Tx})")

    os.makedirs("road_distance", exists_ok=True)
    os.makedirs("angle_to_road", exists_ok=True)
    for (y, x), (box, geo) in progressbar(tiles.items(), redirect_stdout=True):
        boxb = box.buffer(buffer, join_style=2)
        geob = geo.buffer(buffer, join_style=2)

        print(f"Processing tile {(y, x)}.")
        # print("Downloading road network.")
        network = u.download_roads(geob, epsg)
        # print("Rasterizing road network.")
        burned = u.rasterize_roads(network, boxb, epsg, res, res)

        # print("Computing distance and angle to roads.")
        dist, ang = u.dist_ang(burned, res, res)

        profile = dict(
            crs=pyproj.CRS.from_epsg(epsg),
            transform=rasterio.transform.from_bounds(*box.bounds, nx, ny),
        )
        u.save_geotiff(f"road_distance/x{x}_y{y}.tif", dist[mask], **profile)
        u.save_geotiff(f"angle_to_road/x{x}_y{y}.tif", ang[mask], **profile)
