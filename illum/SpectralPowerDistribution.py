#!/usr/bin/env python3

import datetime
import os
from dataclasses import dataclass

import matplotlib as mpl
import numpy as np
import xmltodict

_units_table = {
    "flux": "W / nm",
    "absorptance": "",
    "transmittance": "",
    "reflectance": "",
    "intensity": "W / sr nm",
    "irradiance": "W / m^2 nm",
    "radiance": "W / m^2 nm",
    "exitance": "W / m^2 nm",
    "R-factor": "1 / nm",
    "T-factor": "1 / nm",
    "relative": "",
    "other": "",
}


@dataclass
class SpectralPowerDistribution:
    wavelengths: np.ndarray
    data: np.ndarray
    description: str
    quantity: str

    def __repr__(self):
        return (
            f"SPD<{self.description} "
            f"[{self.wavelengths.min()}:{self.wavelengths.max()}] "
            f"{self.units()}>"
        )

    def units(self):
        return _units_table[self.quantity]

    def interpolate(self, *args, **kwargs):
        return interpolate(self, *args, **kwargs)

    def normalize(self, *args, **kwargs):
        return normalize(self, *args, **kwargs)

    def plot(self, *args, **kwargs):
        return plot(self, *args, **kwargs)

    def to_txt(self, filename, *args, **kwargs):
        return to_txt(filename, self, *args, **kwargs)

    def to_spdx(self, filename, *args, **kwargs):
        return to_spdx(filename, self, *args, **kwargs)


def from_spdx(filename):
    with open(filename) as f:
        content = xmltodict.parse(f.read())["IESTM2714"]

    description = content["Header"]["Description"]
    quantity = content["SpectralDistribution"]["SpectralQuantity"]
    data = [
        [sd["@wavelength"], sd["#text"]]
        for sd in content["SpectralDistribution"]["SpectralData"]
    ]
    data_fp = np.array(data, dtype="float32")

    return SpectralPowerDistribution(
        wavelengths=data_fp[:, 0],
        data=data_fp[:, 1],
        description=description,
        quantity=quantity,
    )


def to_spdx(filename, spd):
    if os.path.splitext(filename)[1].lower() != ".spdx":
        filename += ".spdx"

    content = dict()

    root = content["IESTM2714"] = dict()
    root["@xmlns"] = "iestm2714"
    root["@version"] = "1.0"

    header = root["Header"] = dict()
    header["Description"] = spd.description
    header["DocumentCreator"] = "Illumina"
    timestamp = datetime.datetime.utcnow().isoformat(timespec="seconds")
    header["DocumentCreationDate"] = timestamp

    spct_distr = root["SpectralDistribution"] = dict()
    spct_distr["SpectralQuantity"] = spd.quantity
    spct_distr["SpectralData"] = [
        {"@wavelength": str(wl), "#text": str(val)}
        for wl, val in zip(spd.wavelengths, spd.data)
    ]

    with open(filename, "w") as f:
        xmltodict.unparse(content, f, pretty=True)


def from_txt(filename, skiprows=1, **kwargs):
    data = np.loadtxt(filename, skiprows=skiprows, **kwargs)

    return SpectralPowerDistribution(
        wavelengths=data[:, 0],
        data=data[:, 1],
        description=os.path.basename(filename),
        quantity="other",
    )


def to_txt(filename, spd, /, *, sep="  ", header=True):
    if not os.path.splitext(filename)[1]:
        filename += ".spct"

    if header is True:
        header = ["wavelength", "relativeIntensity"]

    with open(filename, "w") as f:
        if len(header) == 2:
            f.write(sep.join(header) + "\n")
        for wl, val in zip(spd.wavelengths, spd.data):
            f.write(sep.join([str(wl), str(val)]) + "\n")


def from_aster(filename):
    data = np.loadtxt(filename)

    return SpectralPowerDistribution(
        wavelengths=data[:, 0] * 1000,
        data=data[:, 1] / 100,
        description=os.path.basename(filename),
        quantity="reflectance",
    )


def from_file(filename, /):
    name, ext = os.path.splitext(filename)
    funcs = dict(spdx=from_spdx, spct=from_txt, aster=from_aster)
    try:
        return funcs[ext[1:].lower()](filename)
    except KeyError:
        raise ValueError(f"Unknown file type '{ext}'.")


def interpolate(spd, wavelengths):
    # Based on work taken from http://indico.hep.manchester.ac.uk/getFile.py/access?resId=0&materialId=slides&confId=4586
    if type(wavelengths) == SpectralPowerDistribution:
        wavelengths = wavelengths.wavelengths

    def diff(arr):
        return np.diff(arr, prepend=2 * arr[0] - arr[1])

    data = diff(
        np.interp(
            wavelengths,
            spd.wavelengths,
            np.cumsum(diff(spd.wavelengths) * spd.data),
            left=0,
            right=0,
        )
    ) / diff(wavelengths)

    return SpectralPowerDistribution(
        wavelengths=wavelengths,
        data=data,
        description=spd.description,
        quantity=spd.quantity,
    )


def integral(spd, norm=None):
    if norm is None:
        return np.trapz(spd.data, spd.wavelengths)
    return np.trapz(spd.interpolate(norm).data * norm.data, norm.wavelengths)


def normalize(spd, norm=None):
    if norm is None:
        data = spd.data / np.max(spd.data)
        quantity = "relative"
    else:
        data = spd.data / integral(spd, norm=norm)
        quantity = "other"

    return SpectralPowerDistribution(
        wavelengths=spd.wavelengths,
        data=data,
        description=spd.description,
        quantity=quantity,
    )


def plot(spd, /, ax=None, *, axis_labels=False, **kwargs):
    if ax is None:
        ax = mpl.pyplot.gca()

    if axis_labels:
        ax.set_xlabel("Wavelength [nm]")
        ax.set_ylabel(
            f"{spd.quantity}".capitalize()
            + ("" if not spd.units() else f" [{spd.units()}]")
        )
    return ax.plot(spd.wavelengths, spd.data, **kwargs)
