#!/usr/bin/env python3
#
# Inventory preprocessing for Illumina
#
# Author : Alexandre Simoneau
#
# December 2021

import os

import numpy as np
import pandas as pd

import illum


def input_pts(inv_name, params_file="inputs_params.in", dir_name="Inputs"):
    print("Building inputs from discrete inventory.")

    inputs = illum.utils.inputs.prep_inputs(params_file, dir_name)

    inv = pd.read_csv(inv_name, delimiter=" ", dtype={"spct": str, "lop": str})
    domain = illum.PA.load("domain.parr")

    print("Classifying points.")

    coords = illum.utils.geo.transform(t_crs=domain.crs)(inv.lon, inv.lat)
    coords = illum.utils.coords.cart2pol(coords[::-1], domain.center)
    coords = illum.utils.coords.pol2idx(coords, domain.shape, domain._lims)
    inv["y"], inv["x"] = np.array(coords, dtype=int)

    valid = (inv.x >= 0) & (inv.x < domain.shape[1])
    inv = inv[valid].reset_index()

    idx = inv.groupby(["y", "x"]).indices

    print("Calculating the generalized lamps.")

    sources = pd.unique(inv["lop"])

    for s in sources:
        inputs["lops"][s].to_txt(os.path.join(dir_name, f"fctem_{s}"))

    with open(os.path.join(dir_name, "lamps.lst"), "w") as zfile:
        zfile.write("\n".join(sources) + "\n")

    geometry = dict()
    for geo in ["hobs", "dobs", "fobs", "hlmp", "lights"]:
        geometry[geo] = illum.PA.load("domain.parr")

    lumlp = dict()
    for s in sources:
        for wav in inputs["wl"]:
            lumlp[s, wav] = illum.PA.load("domain.parr")

    for (y, x), ind in idx.items():
        lumens = inv["pow"][ind]

        for geo in ["hobs", "dobs", "fobs", "hlmp"]:
            geometry[geo].data[y, x] = np.average(inv[geo][ind], weights=lumens)
        geometry["lights"].data[y, x] = 1

        for s in pd.unique(inv["lop"][ind]):
            mask = inv["lop"][ind] == s
            fctem = np.array(
                [inputs["spcts"][type].data for type in inv["spct"][ind][mask]]
            )
            fctem = np.sum(fctem * lumens[mask].to_numpy()[:, None], 0)

            mean = [np.mean(fctem[mask]) for mask in inputs["bins"]]

            for i, wav in enumerate(inputs["wl"]):
                lumlp[s, wav].data[y, x] = mean[i]

    print("Saving data.")

    for geo, ds in geometry.items():
        ds.save(os.path.join(dir_name, geo))

    for (s, wav), ds in lumlp.items():
        ds.save(os.path.join(dir_name, "%03d_lumlp_%s" % (wav, s)))
