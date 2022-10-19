#!/usr/bin/env python3
#
# Preprocessing for Illumina
#
# Author : Alexandre Simoneau
# modified by Hector Linares for AOD definition
#
# March 2021

import os
import shutil
from glob import glob
from importlib.resources import path

import numpy as np
import yaml
from scipy.interpolate import griddata

import illum
import illum.AngularPowerDistribution as APD
import illum.SpectralPowerDistribution as SPD


def inputs():
    print("Preparing the inputs for the experiment.")

    dir_name = "Inputs/"
    shutil.rmtree(dir_name, True)
    os.makedirs(dir_name)
    shutil.copy("inputs_params.in", dir_name + "inputs_params.in")

    with open("inputs_params.in") as f:
        params = yaml.safe_load(f)

    if params["road_orientation"]:
        print("Computing road orientation")

        with open("domain.in") as f:
            domain_params = yaml.safe_load(f)

        lat = domain_params["latitude"]
        lon = domain_params["longitude"]

        dst, brg = illum.utils.nearest_road(lat, lon)
        np.savetxt(dir_name + "/brng.lst", brg, fmt="%g")

    print("Loading photometry files.")

    # Angular distribution (normalised to 1)
    def parse_key(fname):
        return os.path.basename(fname).rsplit(".", 1)[0].split("_", 1)[0]

    angles = np.arange(181)
    lop = {
        parse_key(fname): APD.from_txt(fname).normalize()
        for fname in glob("Lights/*.lop") + glob("Lights/*.LOP")
    }
    lop.update(
        {
            parse_key(fname): APD.from_ies(fname).normalize()
            for fname in glob("Lights/*.ies") + glob("Lights/*.IES")
        }
    )
    lop = {key: apd.interpolate(step=1) for key, apd in lop.items()}

    # Spectral distribution (normalised with scotopric vision to 1)
    norm_spectrum = SPD.from_txt("Lights/photopic.dat").normalize()
    wav = norm_spectrum.wavelengths
    viirs = (
        SPD.from_txt("Lights/viirs.dat").interpolate(norm_spectrum).normalize()
    )

    spct = {
        parse_key(fname): SPD.from_txt(fname)
        for fname in glob("Lights/*.spct") + glob("Lights/*.SPCT")
    }
    spct.update(
        {
            parse_key(fname): SPD.from_spdx(fname)
            for fname in glob("Lights/*.spdx") + glob("Lights/*.SPDX")
        }
    )
    spct = {
        key: spd.interpolate(norm_spectrum).normalize(norm_spectrum)
        for key, spd in spct.items()
    }

    print("Splitting in wavelengths bins.")

    if os.path.isfile("spectral_bands.dat"):
        bins = np.loadtxt("spectral_bands.dat", delimiter=",")
        n_bins = bins.shape[0]
    else:
        n_bins = params["nb_bins"]
        lmin = params["lambda_min"]
        lmax = params["lambda_max"]

        limits = np.linspace(lmin, lmax, n_bins + 1)
        bins = np.stack([limits[:-1], limits[1:]], axis=1)

    bool_array = (wav >= bins[:, 0:1]) & (wav < bins[:, 1:2])
    x = bins.mean(1)
    bw = bins[:, 1] - bins[:, 0]

    print("Interpolating reflectance.")

    aster = {
        parse_key(fname): SPD.from_aster(fname).interpolate(wav)
        for fname in glob("Lights/*.aster") + glob("Lights/*.ASTER")
    }

    sum_coeffs = sum(
        params["reflectance"][type] for type in params["reflectance"]
    )
    if sum_coeffs == 0:
        sum_coeffs = 1.0

    refl = sum(
        aster[type].data * coeff / sum_coeffs
        for type, coeff in params["reflectance"].items()
    )

    reflect = [np.mean(refl[mask]) for mask in bool_array]

    with open(dir_name + "/refl.lst", "w") as zfile:
        zfile.write("\n".join(["%.06g" % n for n in reflect]) + "\n")

    print("Linking mie files.")

    illumpath = path("illum", "data").as_posix()
    shutil.copy2(
        os.path.abspath(illumpath + "/Molecular_optics/MolecularAbs.txt"),
        dir_name,
    )

    illum.utils.Oillum.PolarArrayC(x)

    shutil.copy("topo.parr", dir_name)

    with open(dir_name + "/wav.lst", "w") as zfile:
        zfile.write("".join("%g %g\n" % (w, b) for w, b in zip(x, bw)))

    if params["zones_inventory"] is not None:
        inv_name = params["zones_inventory"]
        n_inv = 7
        illum.utils.from_zones(
            dir_name,
            inv_name,
            n_inv,
            n_bins,
            params,
            x,
            lop,
            angles,
            wav,
            spct,
            viirs,
            refl,
            bool_array,
        )

    if params["lamps_inventory"] is not None:
        illum.utils.from_lamps(
            dir_name,
            n_bins,
            params,
            x,
            lop,
            angles,
            wav,
            spct,
            viirs,
            refl,
            bool_array,
        )

    # Interpolation of the obstacles properties
    defined = illum.PolarArray.load("domain.parr")
    lights_file = dir_name + "lights.parr"
    if os.path.isfile(lights_file):
        lights = illum.PolarArray.load(lights_file)
        defined.data += lights.data

    for geo in ["obsth", "obstd", "obstf", "altlp"]:
        geometry = illum.PolarArray.load(dir_name + geo + ".parr")
        geometry = (
            griddata(
                points=np.where(defined.data),
                values=geometry[defined.data.astype(bool)],
                xi=tuple(np.ogrid[0 : defined.shape[0], 0 : defined.shape[1]]),
                method="nearest",
            )
            if defined.data.any()
            else np.zeros_like(geometry.data)
        )
        geometry.save(dir_name + geo)

    print("Done.")
