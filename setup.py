#!/usr/bin/env python3

import os
import shutil
from glob import glob

from setuptools import Extension, find_packages, setup
from setuptools.command.build_ext import build_ext


class f2py_Build(build_ext):
    def run(self):
        self.inplace = False
        build_ext.run(self)

    def build_extension(self, ext):
        os.system(
            f"f2py -c {' '.join(ext.sources)} --fcompiler=gfortran"
            f" --f90flags='-fopenmp' -lgomp -m {ext.name}"
        )

        build_py = self.get_finalized_command("build_py")
        src_file, dst_file = self._get_inplace_equivalent(build_py, ext)

        try:
            shutil.copy(src_file, dst_file)
        except (FileNotFoundError, shutil.SameFileError):
            pass


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
