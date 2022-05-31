#!/usr/bin/env python3

import os
import shutil
from glob import glob

import click
import illum


@click.command()
def init():
    """Initialize an execution folder."""
    print("Initializing Illumina execution folder.")
    illumpath = os.path.dirname(illum.__path__[0])

    if not os.path.exists("Lights"):
        shutil.copytree(illumpath + "/Example/Lights", "Lights")

    example_files = glob(illumpath + "/Example/*.*")
    files = [
        s
        for s in example_files
        if not (s.endswith("hdf5") or s.endswith("csv"))
    ]

    for filename in files:
        shutil.copy2(filename, ".")
