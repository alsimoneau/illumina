#!/usr/bin/env python3

import os

import numpy as np

import illum


def spectral_bins(lmin, lmax, N):
    limits = np.linspace(lmin, lmax, N + 1)
    bins = np.stack([limits[:-1], limits[1:]], axis=1)

    return bins


def parse_key(fname):
    return os.path.basename(fname).rsplit(".", 1)[0].split("_", 1)[0]


def open_lops(path):
    return {
        parse_key(fname): illum.APD.from_file(fname)
        .normalize()
        .interpolate(step=1)
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
