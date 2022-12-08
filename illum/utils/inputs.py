#!/usr/bin/env python3

import os as _os
import shutil as _shutil

import numpy as _np
import yaml as _yaml

import illum as _illum


def spectral_bins(lmin, lmax, N):
    limits = _np.linspace(lmin, lmax, N + 1)
    return _np.stack([limits[:-1], limits[1:]], axis=1)


def parse_key(fname):
    return _os.path.basename(fname).rsplit(".", 1)[0].split("_", 1)[0]


def open_lops(path, resolution=1):
    return {
        parse_key(fname): _illum.APD.from_file(fname)
        .normalize()
        .interpolate(step=resolution)
        for fname in _illum.utils.utils.glob_types(
            _os.path.join(path, "*"), ["lop", "ies"]
        )
    }


def open_spcts(path, norm_spct):
    return {
        parse_key(fname): _illum.SPD.from_file(fname)
        .interpolate(norm_spct)
        .normalize(norm_spct)
        for fname in _illum.utils.utils.glob_types(
            _os.path.join(path, "*"), ["spct", "spdx"]
        )
    }


def open_refl(path, norm_spct):
    return {
        parse_key(fname): _illum.SPD.from_aster(fname).interpolate(norm_spct)
        for fname in _illum.utils.utils.glob_types(_os.path.join(path, "*"), ["aster"])
    }


def parse_inventory(filename, n=7):
    """Parse an inventory type file.
    Skips the first 'n' columns."""

    def lamp_norm(lampsData):
        trans = list(map(list, list(zip(*lampsData))))
        norm = sum(trans[0])
        if norm != 0.0:
            trans[0] = [n / norm for n in trans[0]]
        return list(map(list, list(zip(*trans))))

    with open(filename) as inv_file:
        zonData = _illum.utils.utils.strip_comments(inv_file.readlines())
    zonData = [s.split()[n:] for s in zonData]
    zonData = [[s.split("_") for s in i] for i in zonData]
    zonData = [[[float(s[0]), s[1], s[2]] for s in i] for i in zonData]
    zonData = list(map(lamp_norm, zonData))

    return zonData


def prep_inputs(params_file, dir_name):
    _shutil.rmtree(dir_name, True)
    _os.makedirs(dir_name)
    _shutil.copy(params_file, _os.path.join(dir_name, "inputs_params.in"))
    _shutil.copy("topography.parr", _os.path.join(dir_name, "topography.parr"))

    with open(params_file) as f:
        params = _yaml.safe_load(f)

    norm_spct = _illum.SPD.from_txt(_os.path.join("Lights", "photopic.dat"))
    norm_spct = norm_spct.normalize()
    wav = norm_spct.wavelengths
    viirs = _illum.SPD.from_txt(_os.path.join("Lights", "viirs.dat"))
    viirs = viirs.normalize().interpolate(wav)

    lops = _illum.utils.inputs.open_lops("Lights")
    spcts = _illum.utils.inputs.open_spcts("Lights", norm_spct)
    aster = _illum.utils.inputs.open_refl("Lights", norm_spct)

    bins = (
        _np.loadtxt("spectral_bands.dat", delimiter=",")
        if _os.path.isfile("spectral_bands.dat")
        else _illum.utils.inputs.spectral_bins(
            params["lambda_min"], params["lambda_max"], params["nb_bins"]
        )
    )

    bool_array = (wav >= bins[:, 0:1]) & (wav < bins[:, 1:2])
    wl = bins.mean(1)

    sum_coeffs = sum(params["reflectance"][type] for type in params["reflectance"])
    if sum_coeffs == 0:
        sum_coeffs = 1.0

    refl = sum(
        aster[type].data * coeff / sum_coeffs
        for type, coeff in params["reflectance"].items()
    )

    refls = [_np.mean(refl[mask]) for mask in bool_array]

    return dict(lops=lops, spcts=spcts, sens=viirs, wl=wl, bins=bool_array, refls=refls)
