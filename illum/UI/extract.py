#!/usr/bin/env python3
#
# Illumina output extract
#
# Author : Alexandre Simoneau
# unless noted otherwise
#
# February 2019

import os
import re
from collections import defaultdict as ddict
from copy import deepcopy as copy
from functools import partial
from glob import glob

import numpy as np
import pandas as pd

from illum import MultiScaleData as MSD
from illum.pytools import load_bin


def MSDOpen(filename, cached={}):
    if filename in cached:
        return cached[filename]
    ds = MSD.Open(filename)
    cached[filename] = ds
    return ds


def chunker(seq, size):
    return (seq[pos : pos + size] for pos in range(0, len(seq), size))


def add_arrays(a, b):
    if len(a) < len(b):
        a, b = b, a
    c = a.copy()
    c[: len(b)] += b
    return c


def extract(exec_dir, contrib=False, params=(), full=False, profile=False):
    regex_layer = re.compile(r"-layer_(\d+)")
    regex_coords = re.compile(r"observer_coordinates_(-?\d+\.\d+_-?\d+\.\d+)")
    regex_key = re.compile(r"([a-z_]*)_(-?[\d\.]*)")

    parse_key = lambda key: dict(re.findall(regex_key, key))

    skyglow = ddict(float)
    if full:
        outputs = ddict(partial(np.zeros, 6))
    if contrib:
        contributions = dict()
    if profile:
        prof_dist = np.array([])
        profiles = ddict(partial(np.zeros, 0))

    for dirpath, dirnames, filenames in os.walk(exec_dir):
        if not os.path.isfile(os.path.join(dirpath, "illumina.in")):
            continue

        out_names = [fname for fname in filenames if fname.endswith(".out")]

        for oname in out_names:
            if oname == "illumina.out":
                params_name = dirpath.split("exec" + os.sep)[1].replace(os.sep, "-")
            else:
                params_name = oname[:-4]

            key = regex_layer.sub("", params_name)

            try:
                for pname, pvals in params:
                    if pname not in params_name:
                        print("ERROR: Parameter '%s' not found." % pname)
                        exit()
                    for pval in pvals.split(","):
                        if f"{pname}_{pval}" in params_name:
                            break
                    else:
                        raise ValueError()
            except ValueError:
                continue

            with open(os.sep.join([dirpath, oname])) as f:
                lines = f.readlines()
            idx_results = np.where(["==" in line for line in lines])[0][-1] + 1

            val = float(lines[-1])
            if full:
                vals = np.array([float(line) for line in lines[idx_results + 1 :: 2]])
                outputs[key] += vals
            else:
                skyglow[key] += val

            if profile:
                ind = np.where(["==" in line for line in lines])[0]

                prof = np.array(
                    [
                        (
                            float(chunk[2].split()[-2]),
                            float(chunk[3].split()[-2]),
                            float(chunk[-2].split()[-1]),
                        )
                        for chunk in chunker(lines[ind[0] : ind[-1]], ind[1] - ind[0])
                    ]
                )

                profiles[key] = add_arrays(profiles[key], prof[:, 2])

                if len(prof) > len(prof_dist):
                    prof_dist = np.sqrt(prof[:, 0] ** 2 + prof[:, 1] ** 2)

            if contrib:
                try:
                    n_layer = int(regex_layer.search(params_name).groups()[0])
                except AttributeError:
                    # No match, only 1 layer
                    n_layer = 0

                if key not in contributions:
                    try:
                        coords = re.match(regex_coords, params_name).group(1)
                        blank = (
                            dirpath.split("exec")[0]
                            + "/obs_data/%s/blank.hdf5" % coords
                        )
                    except AttributeError:
                        # No match, only 1 coord
                        blank = glob(
                            dirpath.split("exec")[0] + "/obs_data/*/blank.hdf5"
                        )[0]

                    contributions[key] = copy(MSDOpen(blank))

                if oname == "illumina.out":
                    pcl_name = [s for s in filenames if "pcl.bin" in s][0]
                else:
                    pcl_name = "_".join(["pcl", params_name + ".bin"])
                pcl_path = os.path.join(dirpath, pcl_name)
                pcl_data = load_bin(pcl_path)
                pcl_data *= val / pcl_data.sum()
                b = (pcl_data.shape[0] - contributions[key][n_layer].shape[0]) // 2
                contributions[key][n_layer] = pcl_data[b:-b, b:-b] if b else pcl_data

    if full:
        results_names = [s[: s.index("(")].strip() for s in lines[idx_results::2]]
        data = pd.DataFrame.from_dict(
            outputs, "index", columns=results_names
        ).reset_index()
        out = pd.DataFrame(list(data.pop("index").map(parse_key)))
        out = pd.concat([out, data], axis=1).sort_values(list(out.columns))
        out.to_csv("illumina.out", index=False)

        if contrib:
            for key in outputs:
                contributions[key].save(key)
    else:
        data = pd.DataFrame.from_dict(skyglow, "index").reset_index()
        out = pd.DataFrame(list(data.pop("index").map(parse_key)))
        out["Total Radiance"] = data[0]
        out = out.sort_values(list(out.columns))
        out.to_csv("illumina.out", index=False)

        if contrib:
            for key in skyglow:
                contributions[key].save(key)

    if profile:
        length = np.diff(prof_dist, prepend=0)
        for key, val in profiles.items():
            val /= length[: len(val)]
            out = np.stack([prof_dist[: len(val)], val]).T
            np.savetxt(f"profile_{key}.dat", out)
