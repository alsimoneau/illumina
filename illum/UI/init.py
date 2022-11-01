#!/usr/bin/env python3

import os
import shutil
from glob import glob

import illum


def init():
    if illum.path.rsplit(os.sep, 1)[0] in os.path.abspath(os.curdir):
        print("In ILLUMINA code folder. Aborting.")
        return

    if os.path.exists("Lights"):
        print("Existing data files detected. Aborting.")
        return

    print("Initializing ILLUMINA execution folder.")

    shutil.copytree(
        os.path.join(illum.path, "data", "Example", "Lights"), "Lights"
    )

    example_files = glob(os.path.join(illum.path, "data", "Example", "*.*"))
    files = [
        s
        for s in example_files
        if not (s.endswith("hdf5") or s.endswith("csv"))
    ]

    for filename in files:
        shutil.copy2(filename, ".")
