#!/usr/bin/env python3

import warnings

import numpy as np
import rasterio as rio
import rasterio.transform
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
    Na = domain["Nangles"]

    epsg = u.estimate_utm_epsg(lon, lat)

    Nr = round(Na * np.log(rmax) / 2 / np.pi)
    x, y = u.transform(t_crs=epsg)(lon, lat)
    transform = rio.transform.from_origin(
        x - rmax, y + rmax, 2 * rmax, 2 * rmax
    )

    arr = np.zeros((1, 1))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        parr = PA.from_array(
            arr,
            (Na, Nr),
            center=(y, x),
            rmin=rmin,
            rmax=rmax,
            crs=epsg,
            transform=transform,
        )

    parr.save("domain")
