#!/usr/bin/env python3
#
# Inventory preprocessing for Illumina
#
# Author : Alexandre Simoneau
#
# December 2021

import numpy as np

import illum.MultiScaleData as MSD


def input_pts(
    dir_name,
    n_bins,
    params,
    out_name,
    x,
    lop,
    angles,
    wav,
    spct,
    viirs,
    refl,
    bool_array,
):
    print("Building inputs from discrete inventory.")

    # lamps distribution
    inv_name = params["lamps_inventory"]
    lampsData = np.loadtxt(inv_name, usecols=list(range(7)), ndmin=2)
    photometry = np.loadtxt(inv_name, usecols=[-2, -1], dtype=str, ndmin=2)
    domain = MSD.from_domain("domain.ini")

    sources = np.unique(photometry[:, 1])

    print("Classifying points.")

    points = dict()
    for layer in range(len(domain)):
        ysize, xsize = domain[layer].shape
        col, row = domain._get_col_row(lampsData[:, :2].T, layer)
        valid = (0 <= col) * (col < xsize) * (0 <= row) * (row < ysize) > 0
        ind = np.where(valid)[0]
        points[layer] = (col[ind], row[ind], ind)

    print("Calculating the generalized lamps.")

    for n in range(n_bins):
        for s in sources:
            profile = lop[s].vertical_profile()[::-1]
            np.savetxt(
                dir_name + "fctem_wl_%g_lamp_%s.dat" % (x[n], s),
                np.concatenate([profile, angles]).reshape((2, -1)).T,
            )

    with open(dir_name + "lamps.lst", "w") as zfile:
        zfile.write("\n".join(sources) + "\n")

    geometry = dict()
    for geo in ["obsth", "obstd", "obstf", "altlp", "lights"]:
        geometry[geo] = MSD.from_domain("domain.ini")

    lumlp = dict()
    for s in sources:
        for wl in x:
            lumlp[s, wl] = MSD.from_domain("domain.ini")

    for layer, pts in points.items():
        cols, rows, inds = pts
        if len(inds):
            for col, row in np.unique([cols, rows], axis=1).T:
                ind = inds[np.logical_and(cols == col, rows == row)]
                lumens = lampsData[:, 2][ind]

                for n, geo in zip(
                    range(3, 7), ["obsth", "obstd", "obstf", "altlp"]
                ):
                    geometry[geo][layer][row, col] = np.average(
                        lampsData[:, n][ind], weights=lumens
                    )
                geometry["lights"][layer][row, col] = 1

                local_sources = np.unique(photometry[ind][:, 1])
                for s in local_sources:
                    mask = photometry[:, 1][ind] == s
                    fctem = np.array(
                        [
                            spct[type].data
                            for type in photometry[:, 0][ind][mask]
                        ]
                    )
                    fctem = np.sum(fctem * lumens[mask, None], 0)

                    y = [np.mean(fctem[mask]) for mask in bool_array]

                    for i, wl in enumerate(x):
                        lumlp[s, wl][layer][row, col] = y[i]

    print("Saving data.")

    for geo, ds in geometry.items():
        ds.save(dir_name + out_name + "_" + geo)

    for key, ds in lumlp.items():
        s, wl = key
        ds.save(dir_name + "%s_%03d_lumlp_%s" % (out_name, wl, s))
