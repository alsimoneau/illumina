#!/usr/bin/env python3

import numpy as np
import yaml

import illum.PolarArray as PA
import illum.utils as u


def domain():
    with open("domain_params.in") as f:
        domain = yaml.safe_load(f)

    obs_lat = domain.pop("latitude")
    obs_lon = domain.pop("longitude")

    try:
        if len(obs_lat) != len(obs_lon):
            print("ERROR: Latitude and Longitude must have the same length.")
            exit()
    except TypeError:  # lat and lon not lists
        obs_lat = [obs_lat]
        obs_lon = [obs_lon]

    epsg = u.estimate_utm_epsg(obs_lon, obs_lat)
    Nr = round(domain["Nangles"] * np.log(domain["maxRadius"]) / 2 / np.pi)

    PA.from_array(
        np.zeros((1, 1)),
        (domain["Nangles"], Nr),
        rmax=domain["maxRadius"],
        crs=epsg,
    ).save("domain")
