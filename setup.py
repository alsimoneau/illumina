#!/usr/bin/env python3

import os
import re
import shutil
from collections import defaultdict
from glob import glob

from setuptools import Extension, find_packages, setup
from setuptools.command.build_ext import build_ext


def sort_dependancies(files):
    def topological_sort(source):
        """perform topo sort on elements.

        :arg source: list of ``(name, set(names of dependancies))`` pairs
        :returns: list of names, with dependancies listed first

        Taken from stackoverflow.com/q/11557241
        """
        pending = [(name, set(deps)) for name, deps in source]
        emitted = []
        while pending:
            next_pending = []
            next_emitted = []
            for entry in pending:
                name, deps = entry
                deps.difference_update(set((name,)), emitted)
                if deps:
                    next_pending.append(entry)
                else:
                    yield name
                    emitted.append(name)
                    next_emitted.append(name)
            if not next_emitted:
                raise ValueError(
                    "cyclic dependancy detected: %s %r" % (name, (next_pending,))
                )
            pending = next_pending
            emitted = next_emitted

    def_pattern = re.compile(r"(?:FUNCTION|MODULE|SUBROUTINE) ([A-Z_0-9]+)")
    useFct = re.compile(r"([A-Z_0-9]+)\(")
    useMod = re.compile(r"USE ([A-Z_0-9]+)")
    useSub = re.compile(r"CALL ([A-Z_0-9]+)\(")

    def_names = dict()
    use_names = defaultdict(set)

    for fname in files:
        with open(fname) as f:
            content = f.read().upper()

        for match in re.findall(def_pattern, content):
            def_names[match] = fname

        for pattern in (useFct, useSub, useMod):
            for match in re.findall(pattern, content):
                use_names[match].add(fname)

    use_names = {key: val for key, val in use_names.items() if key in def_names}

    depends = {fname: set() for fname in files}
    for name, fnames in use_names.items():
        for fname in fnames:
            depends[fname].add(name)

    depends = {
        fname: {def_names[name] for name in names if def_names[name] != fname}
        for fname, names in depends.items()
    }

    ordered = list(topological_sort(depends.items()))

    return reversed(ordered)


class f2py_Build(build_ext):
    def run(self):
        self.inplace = False
        build_ext.run(self)

    def build_extension(self, ext):
        os.system(
            "f2py -c"
            f" {' '.join(sort_dependancies(ext.sources))} --fcompiler=gfortran"
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
        "joblib",
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
    ext_modules=[Extension("illum.compute.compute", glob("illum/compute/*.f90"))],
    cmdclass=dict(build_ext=f2py_Build),
    include_package_data=True,
    package_data={"illum": ["compute/*", "data/*"]},
)
