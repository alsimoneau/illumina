#!/usr/bin/env python3

import os
import shutil
import site
from glob import glob

from setuptools import Extension, find_packages, setup
from setuptools.command.build_ext import build_ext


class f2py_Build(build_ext):
    def build_extension(self, ext):
        os.system(f"f2py -c {' '.join(ext.sources)} -m {ext.name}")
        fold = list(filter(os.path.exists, site.getsitepackages()))[0]
        so = glob("illum/compute/*.so")
        for s in so:
            shutil.copy(s, os.path.join(fold, s))


with open("illum/__init__.py") as f:
    info = {}
    for line in f:
        if line.startswith("__version__"):
            exec(line, info)
            break

setup(
    name="illum",
    version=info["__version__"],
    packages=find_packages(),
    install_requires=[
        "astropy",
        "Click",
        "fiona",
        "geopandas",
        "GitPython",
        "h5py",
        "matplotlib",
        "numpy",
        "opencv-python",
        "osmnx",
        "pandas",
        "pillow",
        "progressbar2",
        "pyproj",
        "pyyaml",
        "rasterio",
        "scipy",
        "xmltodict",
    ],
    extras_require={"dev": ["black", "flake8", "ipython", "isort"]},
    entry_points="""
        [console_scripts]
        illum=illum.UI.main:main
    """,
    ext_modules=[
        Extension("illum.compute.compute", glob("illum/compute/*.f90"))
    ],
    cmdclass=dict(build_ext=f2py_Build),
    include_package_data=True,
    package_data={"illum": ["compute/*", "data/*"]},
)
