#!/usr/bin/env python3

import os
import shutil
from glob import glob

import numpy as np
import toml

import illum.AngularPowerDistribution as APD
import illum.MultiScaleData as MSD
import illum.pytools as pt
import illum.SpectralPowerDistribution as SPD
from illum.inventory import from_lamps, from_zones


def alternate(name, zones, lights):
    if zones is None and lights is None:
        print("ERROR: At least one of 'zones' and 'lights' must be provided.")
        raise SystemExit

    dirname = "Inputs_%s/" % name

    if os.path.exists(dirname):
        shutil.rmtree(dirname)
    os.makedirs(dirname)

    with open("inputs_params.toml") as f:
        params = toml.load(f)

    if zones is not None and lights is not None:

        print("Validating the inventories.")

        lamps = np.loadtxt(lights, usecols=[0, 1])
        zones = np.loadtxt(params["inventory"]["zones"], usecols=[0, 1, 2])
        zonData = pt.parse_inventory(zones, 0)

        hasLights = [sum(x[0] for x in z) != 0 for z in zonData]

        circles = MSD.from_domain()
        for dat, b in zip(zones, hasLights):
            circles.set_circle((dat[0], dat[1]), dat[2] * 1000, b)

        zones_ind = MSD.from_domain()
        for i, dat in enumerate(zones, 1):
            zones_ind.set_circle((dat[0], dat[1]), dat[2] * 1000, i)

        failed = set()
        for j, coords in enumerate(lamps, 1):
            for i in range(len(circles)):
                try:
                    col, row = circles._get_col_row(coords, i)
                    if circles[i][row, col] and col >= 0 and row >= 0:
                        zon_ind = zones_ind[i][row, col]
                        failed.add((j, coords[0], coords[1], zon_ind))
                except IndexError:
                    continue

        if len(failed):
            for i, lat, lon, zon_ind in sorted(failed):
                print(
                    "WARNING: Lamp #%d (%.06g,%.06g) falls within non-null zone #%d"
                    % (i, lat, lon, zon_ind)
                )
            raise SystemExit()

    shutil.copy("Inputs/inputs_params.toml", dirname)

    print("\nLoading data...")

    # Angular distribution (normalised to 1)
    def parse_key(fname):
        return os.path.basename(fname).rsplit(".", 1)[0].split("_", 1)[0]

    angles = np.arange(181)
    lop = {
        parse_key(fname): APD.load(fname).normalize().interpolate(step=1)
        for fname in glob("data_files/lop/*")
    }

    # Spectral distribution (normalised with scotopric vision to 1 lm / W)
    norm = SPD.load("data_files/photopic.dat").normalize()
    norm.data *= 683.002
    wav = norm.wavelengths
    viirs = SPD.load("data_files/viirs.dat").interpolate(norm).normalize()

    spct = {
        parse_key(fname): SPD.load(fname).interpolate(norm).normalize(norm)
        for fname in glob("data_files/spct/*")
    }

    # Make bins
    if os.path.isfile("spectral_bands.dat"):
        bins = np.loadtxt("spectral_bands.dat", delimiter=",")
        n_bins = bins.shape[0]
    else:
        n_bins = params["wavelengths"]["nb_bins"]
        lmin = params["wavelengths"]["min"]
        lmax = params["wavelengths"]["max"]

        limits = np.linspace(lmin, lmax, n_bins + 1)
        bins = np.stack([limits[:-1], limits[1:]], axis=1)

    bool_array = (wav >= bins[:, 0:1]) & (wav < bins[:, 1:2])
    x = bins.mean(1).tolist()
    bw = bins[:, 1] - bins[:, 0]

    out_name = params["exp_name"]

    aster = {
        parse_key(fname): SPD.load(fname).interpolate(wav)
        for fname in glob("data_files/refl/*")
    }

    sum_coeffs = sum(params["reflectance"][type] for type in params["reflectance"])
    if sum_coeffs == 0:
        sum_coeffs = 1.0

    refl = sum(
        aster[type].data * coeff / sum_coeffs
        for type, coeff in params["reflectance"].items()
    )

    reflect = [np.mean(refl[mask]) for mask in bool_array]
    nspct = [np.mean(norm.data[mask] for mask in bool_array)]

    with open(dirname + "/refl.lst", "w") as zfile:
        zfile.write("\n".join(["%.06g" % n for n in reflect]) + "\n")

    for aero_file in glob("Inputs/*.txt"):
        shutil.copy(aero_file, aero_file.replace("Inputs", dirname))

    shutil.copy("srtm.hdf5", dirname)

    with open(dirname + "/wav.lst", "w") as zfile:
        zfile.write("".join(f"{w:g} {b:g}\n" for w, b in zip(x, bw)))

    if params["inventory"]["zones"]:
        dir_name = ".Inputs_zones/"
        inv_name = params["inventory"]["zones"]
        n_inv = 7
        shutil.rmtree(dir_name, True)
        os.makedirs(dir_name)
        from_zones(
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
        )

        oldlumlp = MSD.from_domain()
        for fname in glob("Inputs/*lumlp*"):
            ds = MSD.Open(fname)
            wl = int(fname.split("_")[1])
            for i, dat in enumerate(ds):
                oldlumlp[i] += dat * nspct[x.index(wl)]

        newlumlp = MSD.from_domain()
        for fname in glob(os.path.join(dir_name, "*lumlp*")):
            ds = MSD.Open(fname)
            wl = int(fname.split("_")[2])
            for i, dat in enumerate(ds):
                newlumlp[i] += dat * nspct[x.index(wl)]

        ratio = MSD.from_domain()
        for i in range(len(ratio)):
            ratio[i] = pt.safe_divide(oldlumlp[i], newlumlp[i])

        for fname in glob(os.path.join(dir_name, "*lumlp*")):
            ds = MSD.Open(fname)
            for i, dat in enumerate(ratio):
                ds[i] *= dat
            ds.save(fname)

    if params["inventory"]["lamps"]:
        dir_name = ".Inputs_lamps/"
        shutil.rmtree(dir_name, True)
        os.makedirs(dir_name)
        from_lamps(
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
        )

    print("Unifying inputs.")

    lfiles = {fname.split(os.sep)[-1] for fname in glob(".Inputs_lamps/*")}
    zfiles = {fname.split(os.sep)[-1] for fname in glob(".Inputs_zones/*")}
    for fname in lfiles - zfiles:
        shutil.move(os.path.join(".Inputs_lamps", fname), dirname)
    for fname in zfiles - lfiles:
        shutil.move(os.path.join(".Inputs_zones", fname), dirname)
    for fname in zfiles & lfiles:
        if "fctem" in fname:
            shutil.move(os.path.join(".Inputs_lamps", fname), dirname)
        elif fname.endswith(".lst"):
            with open(os.path.join(".Inputs_lamps", fname)) as f:
                ldat = f.readlines()
            with open(os.path.join(".Inputs_zones", fname)) as f:
                zdat = f.readlines()
            with open(os.path.join(dirname, fname), "w") as f:
                f.write("".join(sorted(set(ldat + zdat))))
        elif fname.endswith(".hdf5"):
            ldat = MSD.Open(os.path.join(".Inputs_lamps", fname))
            zdat = MSD.Open(os.path.join(".Inputs_zones", fname))
            for i, dat in enumerate(ldat):
                zdat[i][dat != 0] = dat[dat != 0]
            zdat.save(os.path.join(dirname, fname))
        else:
            print("WARNING: File %s not merged properly." % fname)
    if "origin.hdf5" not in zfiles:
        origin = MSD.from_domain()
        origin.save(dirname + "/origin")
    shutil.rmtree(".Inputs_lamps", True)
    shutil.rmtree(".Inputs_zones", True)

    print("Done.")
