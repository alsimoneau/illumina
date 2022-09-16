#!/usr/bin/env python3
#
# Inventory preprocessing for Illumina
#
# Author : Alexandre Simoneau
#
# December 2021

import numpy as np

import illum.compute
import illum.MultiScaleData as MSD
import illum.utils as pt


def from_lamps(
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


def from_zones(
    dir_name,
    inv_name,
    n_inv,
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
    print("Building inputs from zones inventory.")

    # lamps distribution
    zonData = pt.parse_inventory(inv_name, n_inv)

    sources = np.unique([lamp[2] for zd in zonData for lamp in zd])

    for n in range(n_bins):
        for s in sources:
            profile = lop[s].vertical_profile()[::-1]
            np.savetxt(
                dir_name + "fctem_wl_%g_lamp_%s.dat" % (x[n], s),
                np.concatenate([profile, angles]).reshape((2, -1)).T,
            )

    with open(dir_name + "lamps.lst", "w") as zfile:
        zfile.write("\n".join(sources) + "\n")

    print("Making zone property files.")

    circles = MSD.from_domain("domain.ini")  # Same geolocalisation
    zonfile = np.loadtxt(
        params["zones_inventory"], usecols=list(range(7)), ndmin=2
    )

    # zone number
    for i, dat in enumerate(zonfile, 1):
        circles.set_circle((dat[0], dat[1]), dat[2] * 1000, i)
    circles.save(dir_name + out_name + "_zone")

    weights = [sum(z[0] for z in zone) for zone in zonData]
    for w, dat in zip(weights, zonfile):
        circles.set_circle((dat[0], dat[1]), dat[2] * 1000, bool(w))
    circles.save(dir_name + "origin")

    for n, name in zip(range(3, 7), ["obsth", "obstd", "obstf", "altlp"]):
        for i, dat in enumerate(zonfile, 1):
            circles.set_circle((dat[0], dat[1]), dat[2] * 1000, dat[n])
        circles.save(dir_name + out_name + "_" + name)

    print("Inverting lamp intensity.")

    viirs_dat = MSD.Open("stable_lights.hdf5")
    for i in range(len(viirs_dat)):
        viirs_dat[i] *= 1e-5  # nW/cm^2/sr -> W/m^2/sr
        viirs_dat[i][viirs_dat[i] < 0] = 0.0

    water_mask = MSD.Open("water_mask.hdf5")
    for i, wm in enumerate(water_mask):
        viirs_dat[i][wm == 0] = 0.0

    circles = MSD.Open(dir_name + out_name + "_zone.hdf5")
    zon_mask = np.empty(len(circles), dtype=object)
    for i in range(len(zon_mask)):
        zon_mask[i] = (
            np.arange(1, len(zonfile) + 1)[:, None, None] == circles[i]
        )

    spct_keys = list(spct.keys())
    lop_keys = list(lop.keys())

    zones_inventory = np.array(
        [
            [i, lamp[0], spct_keys.index(lamp[1]), lop_keys.index(lamp[2])]
            for i, zon in enumerate(zonData)
            for lamp in zon
        ]
    )
    zones_inventory[:, [0, 2, 3]] += 1  # Fortran indexing

    spcts = np.array(list(spct.values()))
    lops = np.array(list(lop.values()))
    sources_idx = [
        list(sources).index(k) if k in sources else -1 for k in lop_keys
    ]
    sources_idx = np.array(sources_idx) + 1

    lumlps = [
        illum.compute.viirs2lum(
            nzones=len(zonData),
            nsources=len(sources),
            viirs=viirs_dat[i],
            zones=circles[i],
            angles=angles,
            wav=wav,
            bands=bool_array,
            sens=viirs,
            lops=np.array([lops[k].vertical_profile() for k in lop_keys]),
            spcts=np.array([spcts[k].data for k in spct_keys]),
            sources=sources_idx,
            ivtr=zones_inventory,
            pixsize=circles.pixel_size(i),
            reflect=refl,
        )
        for i in range(len(circles))
    ]

    for n, wl in enumerate(x):
        for s, source in enumerate(sources):
            new = MSD.from_domain("domain.ini")
            for layer in range(len(circles)):
                new[layer] = lumlps[layer][n, s]
            new.save(dir_name + "%s_%g_lumlp_%s" % (out_name, wl, source))
