#!/usr/bin/env python3

import os
import shutil

import numpy as np
import yaml

import illum


def spectral_bins(lmin, lmax, N):
    limits = np.linspace(lmin, lmax, N + 1)
    return np.stack([limits[:-1], limits[1:]], axis=1)


def parse_key(fname):
    return os.path.basename(fname).rsplit(".", 1)[0].split("_", 1)[0]


def open_lops(path, resolution=1):
    return {
        parse_key(fname): illum.APD.from_file(fname)
        .normalize()
        .interpolate(step=resolution)
        for fname in illum.utils.glob_types(
            os.path.join(path, "*"), ["lop", "ies"]
        )
    }


def open_spcts(path, norm_spct):
    return {
        parse_key(fname): illum.SPD.from_file(fname)
        .interpolate(norm_spct)
        .normalize(norm_spct)
        for fname in illum.utils.glob_types(
            os.path.join(path, "*"), ["spct", "spdx"]
        )
    }


def open_refl(path, norm_spct):
    return {
        parse_key(fname): illum.SPD.from_aster(fname).interpolate(norm_spct)
        for fname in illum.utils.glob_types(os.path.join(path, "*"), ["aster"])
    }


def parse_inventory(filename, n=0):
    """Parse an inventory type file.
    Skips the first 'n' columns."""

    def lamp_norm(lampsData):
        trans = list(map(list, list(zip(*lampsData))))
        norm = sum(trans[0])
        if norm != 0.0:
            trans[0] = [n / norm for n in trans[0]]
        return list(map(list, list(zip(*trans))))

    with open(filename) as inv_file:
        zonData = illum.utils.strip_comments(inv_file.readlines())
    zonData = [s.split()[n:] for s in zonData]
    zonData = [[s.split("_") for s in i] for i in zonData]
    zonData = [[[float(s[0]), s[1], s[2]] for s in i] for i in zonData]
    zonData = list(map(lamp_norm, zonData))

    return zonData


def prep_inputs(params_file, dir_name):
    shutil.rmtree(dir_name, True)
    os.makedirs(dir_name)
    shutil.copy(params_file, os.path.join(dir_name, "inputs_params.in"))

    with open(params_file) as f:
        params = yaml.safe_load(f)

    norm_spct = illum.SPD.from_txt(os.path.join("Lights", "photopic.dat"))
    norm_spct = norm_spct.normalize()
    wav = norm_spct.wavelengths
    viirs = illum.SPD.from_txt(os.path.join("Lights", "viirs.dat"))
    viirs = viirs.normalize().interpolate(wav)

    lops = illum.utils.open_lops("Lights")
    spcts = illum.utils.open_spcts("Lights", norm_spct)
    aster = illum.utils.open_refl("Lights", norm_spct)

    bins = (
        np.loadtxt("spectral_bands.dat", delimiter=",")
        if os.path.isfile("spectral_bands.dat")
        else illum.utils.spectral_bins(
            params["lambda_min"], params["lambda_max"], params["nb_bins"]
        )
    )

    bool_array = (wav >= bins[:, 0:1]) & (wav < bins[:, 1:2])
    wl = bins.mean(1)

    sum_coeffs = sum(
        params["reflectance"][type] for type in params["reflectance"]
    )
    if sum_coeffs == 0:
        sum_coeffs = 1.0

    refl = sum(
        aster[type].data * coeff / sum_coeffs
        for type, coeff in params["reflectance"].items()
    )

    refls = [np.mean(refl[mask]) for mask in bool_array]

    return lops, spcts, viirs, wl, bool_array, refls
