#!/usr/bin/env python3

import importlib.resources
import os
import shutil
from glob import glob


def init():
    illumpath = importlib.resources.files("illum").as_posix()
    if illumpath.rsplit(os.sep, 1)[0] in os.path.abspath(os.curdir):
        print("In ILLUMINA code folder. Aborting.")
        return

    if os.path.exists("Lights"):
        print("Existing data files detected. Aborting.")
        return

    print("Initializing ILLUMINA execution folder.")

    shutil.copytree(
        illumpath + os.sep.join(["", "data", "Example", "Lights"]),
        "Lights",
    )

    example_files = glob(
        illumpath + os.sep.join(["", "data", "Example", "*.*"])
    )
    files = [
        s
        for s in example_files
        if not (s.endswith("hdf5") or s.endswith("csv"))
    ]

    for filename in files:
        shutil.copy2(filename, ".")
