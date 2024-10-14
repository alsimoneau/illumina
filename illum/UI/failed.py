#!/usr/bin/env python3

import fnmatch
import os
from glob import glob


def recursive_glob(rootdir=".", pattern="*"):
    for root, dirnames, filenames in os.walk(rootdir):
        for filename in fnmatch.filter(filenames, pattern):
            yield os.path.join(root, filename)


def failed(scheduler=None):
    if scheduler is not None and scheduler != "none":

        execute = dict(
            parallel="./execute &", sequential="./execute", slurm="sbatch ./execute"
        )
        execute_str = execute[scheduler] + "\n"

        def failed(dirname):
            return "\n".join(
                ["cd " + os.path.abspath(dirname), execute_str, "sleep 0.05"]
            )

    else:

        def failed(dirname):
            return dirname

    out = []
    for dirname in recursive_glob(pattern="illumina"):
        dirname = os.path.dirname(dirname)
        if not os.path.isfile(os.path.join(dirname, "illumina.in")):
            out.append(failed(dirname))
        else:
            with open(os.path.join(dirname, "illumina.in")) as f:
                basename = f.readlines()[1].split()[0]

            outnames = glob(os.path.join(dirname, basename + "*.out"))
            nb_in = len(glob(os.path.join(dirname, "*.in")))

            if nb_in <= 2:
                if len(outnames) == 0:
                    out.append(failed(dirname))
                else:
                    with open(outnames[0]) as f:
                        lines = f.readlines()
                    if len(lines) < 2 or "Diffuse radiance" not in lines[-2]:
                        out.append(failed(dirname))
            elif (
                os.path.isfile(os.path.join(dirname, basename + ".out"))
                or len(outnames) + 1 != nb_in
            ):
                out.append(failed(dirname))

    return out
