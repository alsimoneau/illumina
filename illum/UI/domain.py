#!/usr/bin/env python3

import numpy as np
import rasterio as rio
import yaml

import illum.PolarArray as PA
import illum.utils as u


def domain():
    with open("domain_params.in") as f:
        domain = yaml.safe_load(f)

    lat = domain["latitude"]
    lon = domain["longitude"]
    rmin = domain["minRadius"]
    rmax = domain["maxRadius"]
    res = domain["resolution"]
    Na = domain["Nangles"]

    epsg = u.estimate_utm_epsg(lon, lat)
    x, y = u.transform(s_crs=epsg)(lon, lat)
    n = round(rmax / res) + 0.5
    transform = rio.transform.from_origin(
        west=x - n * res, north=y + n * res, xsize=res, ysize=res
    )

    Nr = round(Na * np.log(rmax) / 2 / np.pi)

    arr = np.zeros((1, 1))
    parr = PA.from_array(
        arr, (Na, Nr), rmin=rmin, rmax=rmax, crs=epsg, transform=transform
    )
    parr.center = (int(n), int(n))
    parr.shape = (int(2 * n), int(2 * n))
    parr.save("domain")
