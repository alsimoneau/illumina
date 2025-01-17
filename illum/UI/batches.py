#!/usr/bin/env python3
#
# batch processing
#
# Author : Alexandre Simoneau
# unless noted otherwise
#
# March 2021

import os
import shutil
from collections import ChainMap, OrderedDict
from functools import partial
from glob import glob
from itertools import product

import numpy as np
import toml
from progressbar import progressbar

import illum
from illum import MultiScaleData as MSD
from illum.pytools import save_bin

progress = partial(progressbar, redirect_stdout=True)


def input_line(val, comment, n_space=30):
    value_str = " ".join(str(v) for v in val)
    comment_str = " ; ".join(comment)
    return "%-*s ! %s" % (n_space, value_str, comment_str)


def MSDOpen(filename, cached={}):
    if filename in cached:
        return cached[filename]
    ds = MSD.Open(filename)
    cached[filename] = ds
    return ds


def batches(
    input_path=".",
    compact=False,
    batch_size=300,
    scheduler="sequential",
    batch_name="batch",
):
    execute = dict(
        parallel="./execute &", sequential="./execute", slurm="sbatch ./execute"
    )
    execute_str = execute[scheduler] + "\n"

    os.chdir(input_path)

    with open("inputs_params.toml") as f:
        params = toml.load(f)

    for fname in glob(batch_name + "*"):
        os.remove(fname)

    ds = MSD.Open(glob("*.hdf5")[0])

    # Pre process the obs extract
    print("Preprocessing...")
    shutil.rmtree("obs_data", True)
    lats, lons = ds.get_obs_pos()
    xs, ys = ds.get_obs_pos(proj=True)
    for lat, lon in zip(lats, lons):
        for i in range(len(ds)):
            os.makedirs("obs_data/%6f_%6f/%d" % (lat, lon, i))

    for i, fname in enumerate(progress(glob("*.hdf5")), 1):
        dataset = MSD.Open(fname)
        for clipped in dataset.split_observers():
            lat, lon = clipped.get_obs_pos()
            lat, lon = lat[0], lon[0]

            if "lumlp" in fname:
                clipped.set_buffer(0)
                clipped.set_overlap(0)
            for i, dat in enumerate(clipped):
                padded_dat = np.pad(dat, (512 - dat.shape[0]) // 2, "constant")
                save_bin(
                    "obs_data/%6f_%6f/%i/%s"
                    % (lat, lon, i, fname.rsplit(".", 1)[0] + ".bin"),
                    padded_dat,
                )
            if "srtm" in fname:
                for j in range(len(clipped)):
                    clipped[j][:] = 0
                clipped.save(f"obs_data/{lat:6f}_{lon:6f}/blank")

    # Add wavelength and multiscale
    spectral_bands = np.loadtxt("wav.lst", ndmin=2)
    params["wavelength"] = spectral_bands[:, 0].tolist()
    params["layer"] = list(range(len(ds)))
    params["observer_coordinates"] = list(zip(*ds.get_obs_pos()))

    wls = params["wavelength"]
    refls = np.loadtxt("refl.lst", ndmin=1).tolist()

    for pname in ["layer", "observer_coordinates"]:
        if len(params[pname]) == 1:
            params[pname] = params[pname][0]

    with open("lamps.lst") as f:
        lamps = f.read().split()

    if os.path.isfile("brng.lst"):
        brng = np.loadtxt("brng.lst", ndmin=1)

    # Clear and create execution folder
    dir_name = "exec" + os.sep
    shutil.rmtree(dir_name, True)
    os.makedirs(dir_name)

    count = 0
    params.update(
        {
            f"{t}_{key}": val
            for t, k in params.items()
            if isinstance(k, dict)
            for key, val in k.items()
        }
    )
    multival = [k for k in params if isinstance(params[k], list) and k != "aerosols"]
    multival = sorted(multival, key=len, reverse=True)  # Semi-arbitrary sort
    param_space = [params[k] for k in multival]

    N = np.prod([len(p) for p in param_space])
    for param_vals in progress(product(*param_space), max_value=N):
        local_params = OrderedDict(zip(multival, param_vals))
        P = ChainMap(local_params, params)
        if (
            "viewing_angles_azimuth" in multival
            and P["viewing_angles_elevation"] == 90
            and params["viewing_angles_azimuth"].index(P["viewing_angles_azimuth"]) != 0
        ):
            continue

        if os.path.isfile("brng.lst"):
            obs_index = (
                0
                if "observer_coordinates" not in multival
                else params["observer_coordinates"].index(P["observer_coordinates"])
            )
            bearing = brng[obs_index]
        else:
            bearing = 0

        coords = "%6f_%6f" % P["observer_coordinates"]
        if "observer_coordinates" in multival:
            P["observer_coordinates"] = coords

        if compact:
            fold_name = (
                dir_name
                + os.sep.join(
                    f"{k}_{v}"
                    for k, v in local_params.items()
                    if k in ["observer_coordinates", "wavelength", "layer"]
                )
                + os.sep
            )
        else:
            fold_name = (
                dir_name
                + os.sep.join(f"{k}_{v}" for k, v in local_params.items())
                + os.sep
            )

        unique_ID = "-".join("%s_%s" % item for item in local_params.items())
        wavelength = "%g" % P["wavelength"]
        layer = P["layer"]
        reflectance = refls[wls.index(P["wavelength"])]
        bandwidth = spectral_bands[wls.index(P["wavelength"]), 1]

        if not os.path.isdir(fold_name):
            os.makedirs(fold_name)
            # Linking files
            for i, aerosols in enumerate(params["aerosols"], 1):
                mie_file = "{}_{}.txt".format(aerosols["profile"], wavelength)
                os.symlink(
                    os.path.relpath(mie_file, fold_name),
                    fold_name + f"aerosols_{i}.txt",
                )

            os.symlink(
                os.path.relpath("MolecularAbs.txt", fold_name),
                fold_name + "MolecularAbs.txt",
            )

            for i, lamp in enumerate(lamps, 1):
                os.symlink(
                    os.path.relpath(
                        f"fctem_wl_{wavelength}_lamp_{lamp}.dat", fold_name
                    ),
                    fold_name + "fctem_%03d.dat" % i,
                )

            illumpath = os.path.dirname(illum.__path__[0])
            os.symlink(
                os.path.abspath(illumpath + "/bin/illumina"), fold_name + "illumina"
            )

            # Copying layer data
            obs_fold = os.path.join("obs_data", coords, str(layer))

            os.symlink(
                os.path.relpath(os.path.join(obs_fold, "srtm.bin"), fold_name),
                fold_name + "topogra.bin",
            )

            os.symlink(
                os.path.relpath(os.path.join(obs_fold, "origin.bin"), fold_name),
                fold_name + "origin.bin",
            )

            for name in ["obstd", "obsth", "obstf", "altlp"]:
                os.symlink(
                    os.path.relpath(os.path.join(obs_fold, f"{name}.bin"), fold_name),
                    fold_name + f"{name}.bin",
                )

            for i, lamp in enumerate(lamps, 1):
                os.symlink(
                    os.path.relpath(
                        os.path.join(obs_fold, f"{wavelength}_lumlp_{lamp}.bin"),
                        fold_name,
                    ),
                    fold_name + f"lumlp_{i:03d}.bin",
                )

        # Create illumina.in
        input_data = [
            (("", "Input file for ILLUMINA"),),
            (("", ""),),
            (
                (ds.pixel_size(layer), "Cell size along X [m]"),
                (ds.pixel_size(layer), "Cell size along Y [m]"),
            ),
            (
                (256, "Observer X position"),
                (256, "Observer Y position"),
                (P["observer_elevation"], "Observer elevation above ground [m]"),
            ),
            ((len(lamps), "Number of source types"),),
            ((wavelength, "Wavelength [nm]"), (bandwidth, "Bandwidth [nm]")),
            (("", ""),),
            (
                (P["viewing_angles_elevation"], "Elevation viewing angle"),
                (
                    (P["viewing_angles_azimuth"] + bearing) % 360,
                    "Azimuthal viewing angle",
                ),
            ),
            ((P["viewing_angles_direct_fov"], "Direct field of view"),),
            (("", ""),),
            ((P["stop_limit"], "Contribution threshold"),),
            ((P["scattering_single"] * 1, "Single scattering activated"),),
            ((P["scattering_double"] * 1, "Double scattering activated"),),
            (("", ""),),
            ((P["observer_obstacles"] * 1, "Obstacles around observer"),),
            ((reflectance, "Reflectance"),),
            (
                (
                    P["reflection_radius"],
                    "Radius around light sources where reflextions are computed",
                ),
            ),
            (
                (
                    P["clouds_model"],
                    "Cloud model: "
                    "0=clear, "
                    "1=Thin Cirrus/Cirrostratus, "
                    "2=Thick Cirrus/Cirrostratus, "
                    "3=Altostratus/Altocumulus, "
                    "4=Cumulus/Cumulonimbus, "
                    "5=Stratocumulus",
                ),
                (P["clouds_base"], "Cloud base altitude [m]"),
                (P["clouds_fraction"], "Cloud fraction"),
            ),
            (("", ""),),
            ((P["air_pressure"], "Ground level pressure [kPa]"),),
            ((len(P["aerosols"]), "Number of aerosols layers"),),
        ]
        for i, aerosols in enumerate(P["aerosols"], 1):
            input_data.append(
                (
                    (f"aerosols_{i}.txt", "Optical cross section file"),
                    (aerosols["optical_depth"], "Optical depth at 500nm"),
                    (aerosols["angstrom_coefficient"], "Angstom coefficient"),
                    (aerosols["scale_heigth"], "Scale height [m]"),
                )
            )

        with open(fold_name + unique_ID + ".in", "w") as f:
            lines = (input_line(*zip(*line_data)) for line_data in input_data)
            f.write("\n".join(lines))

        # Write execute script
        if not os.path.isfile(fold_name + "execute"):
            with open(fold_name + "execute", "w") as f:
                f.write("#!/bin/sh\n")
                f.write("#SBATCH --job-name=Illumina\n")
                f.write("#SBATCH --mem=2G\n")
                f.write("cd %s\n" % os.path.abspath(fold_name))
                f.write("umask 0011\n")
            os.chmod(fold_name + "execute", 0o777)

            # Append execution to batch list
            with open(f"{"batch_name"}_{(count//batch_size)+1}", "a") as f:
                f.write("cd %s\n" % os.path.abspath(fold_name))
                f.write(execute_str)
                f.write("sleep 0.05\n")

            count += 1

        # Add current parameters execution to execution script
        with open(fold_name + "execute", "a") as f:
            f.write("cp %s.in illumina.in\n" % unique_ID)
            f.write("./illumina\n")
            f.write(f"mv illumina.out {unique_ID}.out\n")
            f.write(f"mv illumina_pcl.bin pcl_{unique_ID}.bin\n")

    print("Final count:", count)

    print("Done.")
