#!/usr/bin/env python3

import inspect
import warnings

import numpy as np
import rasterio as rio
import rasterio.transform
import yaml

import illum


def domain(params="domain_params.in"):
    if type(params) == str:
        with open(params) as f:
            domain = yaml.safe_load(f)
    elif type(params) == dict:
        domain = params
    else:
        raise TypeError("'params' must be a filename or a dictionnary.")

    lat = domain["latitude"]
    lon = domain["longitude"]
    rmin = domain["minRadius"]
    rmax = domain["maxRadius"]
    Na = domain["Nangles"]

    epsg = illum.utils.estimate_utm_epsg(lon, lat)

    Nr = round(np.log(rmax / rmin) / np.arcsinh(np.pi / Na) / 2)
    x, y = illum.utils.transform(t_crs=epsg)(lon, lat)
    transform = rio.transform.from_origin(
        x - rmax, y + rmax, 2 * rmax, 2 * rmax
    )

    arr = np.zeros((1, 1))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        parr = illum.PA.from_array(
            arr,
            (Na, Nr),
            center=(y, x),
            rmin=rmin,
            rmax=rmax,
            crs=epsg,
            transform=transform,
        )

    plat = f"{abs(lat):.6f}°{'N' if lat>=0 else 'S'}"
    plon = f"{abs(lon):.6f}°{'E' if lon>=0 else 'W'}"
    prmax = f"{rmax:,}m".replace(",", "'")
    prmin = f"{rmin:,}m".replace(",", "'")
    ra = f"{360/Na:.3g}°"
    rr = f"{rmin*((rmax/rmin)**(1/Nr)-1):.3g}m"

    printout = inspect.cleandoc(
        f"""center coordinate: {plat}, {plon}
        radius limits: {prmin}, {prmax}
        polar resolution: {Na} angles ({ra}), {Nr} radii ({rr})
        crs: EPSG {epsg}"""
    )
    print(printout)

    return parr.save("domain")
