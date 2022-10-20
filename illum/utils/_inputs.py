#!/usr/bin/env python3

import os

import numpy as np

import illum.AngularPowerDistribution as APD
import illum.SpectralPowerDistribution as SPD
import illum.utils as u


def spectral_bins(lmin, lmax, N):
    limits = np.linspace(lmin, lmax, N + 1)
    bins = np.stack([limits[:-1], limits[1:]], axis=1)

    return bins


def parse_key(fname):
    return os.path.basename(fname).rsplit(".", 1)[0].split("_", 1)[0]


def open_lops(path):
    lop = {
        parse_key(fname): APD.from_file(fname).normalize()
        for fname in u.glob_types(os.path.join(path, "*"), ["lop", "ies"])
    }
    return {key: apd.interpolate(step=1) for key, apd in lop.items()}


def open_spcts(path):
    norm_spectrum = SPD.from_txt(
        os.path.join(path, "photopic.dat")
    ).normalize()

    spct = {
        parse_key(fname): SPD.from_file(fname)
        for fname in u.glob_types(os.path.join(path, "*"), ["spct", "spdx"])
    }
    return {
        key: spd.interpolate(norm_spectrum).normalize(norm_spectrum)
        for key, spd in spct.items()
    }
