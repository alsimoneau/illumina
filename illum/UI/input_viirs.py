#!/usr/bin/env python3
#
# Inventory preprocessing for Illumina
#
# Author : Alexandre Simoneau
#
# October 2022

import os

import numpy as np

import illum


def input_viirs(inv_name, params_file="inputs_params.in", dir_name="Inputs"):
    print("Building inputs from zones inventory.")

    inputs = illum.utils.prep_inputs(params_file, dir_name)
    lops, spcts, viirs, wl, bool_array, refls = inputs

    # lamps distribution
    zonData = illum.utils.parse_inventory(inv_name, 7)

    sources = {lamp[2] for zd in zonData for lamp in zd}

    for s in sources:
        lops[s].to_txt(os.path.join(dir_name, f"fctem_{s}"))

    with open(os.path.join(dir_name, "lamps.lst"), "w") as f:
        f.write("\n".join(sources) + "\n")

    print("Making zone property files.")

    circles = illum.PA.load("domain.parr")
    zonfile = np.loadtxt(inv_name, usecols=list(range(7)), ndmin=2)

    # zone number
    for i, dat in enumerate(zonfile, 1):
        circles.set_circle((dat[0], dat[1]), dat[2] * 1000, i)
    circles.save(os.path.join(dir_name, "zone"))

    weights = [sum(z[0] for z in zone) for zone in zonData]
    for w, dat in zip(weights, zonfile):
        circles.set_circle((dat[0], dat[1]), dat[2] * 1000, bool(w))
    circles.save(os.path.join(dir_name, "origin"))

    for n, name in zip(range(3, 7), ["obsth", "obstd", "obstf", "altlp"]):
        for i, dat in enumerate(zonfile, 1):
            circles.set_circle((dat[0], dat[1]), dat[2] * 1000, dat[n])
        circles.save(os.path.join(dir_name, name))

    print("Inverting lamp intensity.")

    viirs_dat = illum.PA.load("lights.parr")
    viirs_dat.data *= 1e-5  # nW/cm^2/sr -> W/m^2/sr
    viirs_dat.data[viirs_dat.data < 0] = 0.0

    water_mask = illum.PA.load("water.parr")
    viirs_dat.data[water_mask.data == 0] = 0.0

    circles = illum.PA.load(os.path.join(dir_name, "zone.parr"))

    spct_keys = list(spcts.keys())
    lop_keys = list(lops.keys())

    zones_inventory = np.array(
        [
            [i, lamp[0], spct_keys.index(lamp[1]), lop_keys.index(lamp[2])]
            for i, zon in enumerate(zonData)
            for lamp in zon
        ]
    )
    zones_inventory[:, [0, 2, 3]] += 1  # Fortran indexing

    sources_idx = [
        list(sources).index(k) if k in sources else -1 for k in lop_keys
    ]
    sources_idx = np.array(sources_idx) + 1

    lumlps = illum.compute.viirs2lum(
        nzones=len(zonData),
        nsources=len(sources),
        viirs=viirs_dat.data,
        zones=circles.data.astype("uint32"),
        angles=lops[lop_keys[0]].vertical_angles,
        wav=spcts[spct_keys[0]].wavelengths,
        bands=bool_array,
        sens=viirs.data,
        lops=np.array([lops[k].vertical_profile() for k in lop_keys]),
        spcts=np.array([spcts[k].data for k in spct_keys]),
        sources=sources_idx,
        ivtr=zones_inventory,
        pixsize=circles.area(),
        reflect=refls,
    )

    for n, x in enumerate(wl):
        for s, source in enumerate(sources):
            new = illum.PA.load("domain.parr")
            new.data = lumlps[n, s]
            new.save(os.path.join(dir_name, "%g_lumlp_%s" % (x, source)))
