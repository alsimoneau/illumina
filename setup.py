#!/usr/bin/env python3

import os

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext


class f2py_Extension(Extension):
    def __init__(self, name):
        Extension.__init__(self, name, sources=[])


class f2py_Build(build_ext):
    def run(self):
        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        name = os.path.basename(ext.name)
        os.system(f"cd {ext.name}; f2py -c  `ls *.f90 | xargs` -m {name}")


with open("illum/__init__.py") as f:
    info = {}
    for line in f:
        if line.startswith("__version__"):
            exec(line, info)
            break

setup(
    name="illum",
    version=info["__version__"],
    packages=["illum"],
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
    ext_modules=[f2py_Extension("illum/compute")],
    cmdclass=dict(build_ext=f2py_Build),
)
