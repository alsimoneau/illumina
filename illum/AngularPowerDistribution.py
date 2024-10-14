#!/usr/bin/env python3

import os
from dataclasses import dataclass

import matplotlib as mpl
import numpy as np
import scipy.interpolate


def mids(arr, /):
    return np.concatenate([[arr[0]], np.mean([arr[1:], arr[:-1]], 0), [arr[-1]]])


@dataclass
class AngularPowerDistribution:
    vertical_angles: np.ndarray
    horizontal_angles: np.ndarray
    data: np.ndarray
    type: int = 1
    lumens: float = -1

    def __repr__(self):
        return f"APD<{self.data.shape}>"

    def cycle(self, *args, **kwargs):
        return cycle(self, *args, **kwargs)

    def interpolate(self, *args, **kwargs):
        return interpolate(self, *args, **kwargs)

    def normalize(self, *args, **kwargs):
        return normalize(self, *args, **kwargs)

    def plotV(self, *args, **kwargs):
        return plotV(self, *args, **kwargs)

    def plotH(self, *args, **kwargs):
        return plotH(self, *args, **kwargs)

    def plot2d(self, *args, **kwargs):
        return plot2d(self, *args, **kwargs)

    def plot3d(self, *args, **kwargs):
        return plot3d(self, *args, **kwargs)

    def save(self, filename, *args, **kwargs):
        return save(filename, self, *args, **kwargs)

    def vertical_profile(self, *args, **kwargs):
        return vertical_profile(self, *args, **kwargs)

    def horizontal_profile(self, *args, **kwargs):
        return horizontal_profile(self, *args, **kwargs)


def load(filename, /, facing=None, **kwargs):
    match os.path.splitext(filename)[1].lower():
        case ".ies":
            apd, f = from_iesna(filename, **kwargs), 0
        case ".cib" | ".tml" | ".c1s" | ".cc":
            apd, f = from_cibse(filename, **kwargs), 90
        case ".ldt" | ".exl":
            apd, f = from_eulumdat(filename, **kwargs), 90
        case _:
            apd, f = from_ascii(filename, **kwargs), 0

    apd.horizontal_angles -= facing if facing is not None else f
    return apd.cycle()


def save(filename, apd, /, **kwargs):
    match os.path.splitext(filename)[1].lower():
        case ".ies":
            to_iesna(filename, **kwargs)
        case ".cib" | ".tml" | ".c1s" | ".cc":
            to_cibse(filename, **kwargs)
        case ".ldt" | ".exl":
            to_eulumdat(filename, **kwargs)
        case _:
            to_ascii(filename, **kwargs)


def from_iesna(filename, /, facing=0):
    with open(filename) as f:
        for line in f:
            if line.startswith("TILT="):
                break
        skip = 0 if "NONE" in line else 4
        data = f.read().replace(",", " ").split()[skip:]
    lumens = int(data[0]) * float(data[1])
    factor = float(data[2])
    nV = int(data[3])
    nH = int(data[4])
    data = data[13:]
    va, data = np.array(data[:nV], dtype="float32"), data[nV:]
    ha, data = np.array(data[:nH], dtype="float32"), data[nH:]
    data = np.reshape(data, (nH, nV)).T.astype("float32") * factor

    return AngularPowerDistribution(
        lumens=lumens, vertical_angles=va, horizontal_angles=ha, data=data
    )


def to_iesna(filename, apd, /, facing=0):
    out = [
        "IES:LM-63-2019",
        "TILT=NONE",
        f"1 {apd.lumens} 1.0 {apd.data.shape[0]} {apd.data.shape[1]} "  # cont.
        f"{apd.type} 1 0 0 0",
        "1.0 1.0 0.0",
    ]
    with open(filename, "w") as f:
        f.write("\n".join(out) + "\n")
        apd.vertical_angles.tofile(f, sep=" ", format="%f")
        f.write("\n")
        apd.horizontal_angles.tofile(f, sep=" ", format="%f")
        f.write("\n")
        for row in apd.data.T:
            row.tofile(f, sep=" ", format="%f")
            f.write("\n")


def from_cibse(filename, /, facing=90):
    with open(filename) as f:
        for i, line in enumerate(f):
            if i == 7:
                break
        data = f.read().split()
        factor = float(data[5])
        nV = int(data[9])
        nH = int(data[10])
        data = data[11:]
        va, data = np.array(data[:nV], dtype="float32"), data[nV:]
        ha, data = np.array(data[:nH], dtype="float32"), data[nH:]
        data = data[: nV * nH]
        data = np.reshape(data, (nH, nV)).T.astype("float32") * factor

    return AngularPowerDistribution(
        lumens=1, vertical_angles=va, horizontal_angles=ha, data=data
    )


def to_cibse(filename, apd, /, facing=90):
    raise NotImplementedError


def from_eulumdat(filename, /, facing=90):
    with open(filename) as f:
        data = f.readlines()
    nH = int(data[3])
    nV = int(data[5])
    factor = float(data[23])
    lumens = float(data[28])
    data = data[42:]
    ha, data = np.array(data[:nH], dtype="float32"), data[nH:]
    va, data = np.array(data[:nV], dtype="float32"), data[nV:]
    data = np.reshape(data[: nH * nV], (nH, nV)).T.astype("float32") * factor

    return AngularPowerDistribution(
        lumens=lumens, vertical_angles=va, horizontal_angles=ha, data=data
    )


def to_eulumdat(filename, apd, /, facing=90):
    raise NotImplementedError


def from_ascii(filename, /):
    data, ang = np.loadtxt(filename)[::-1].T
    return AngularPowerDistribution(
        type=1,
        vertical_angles=180 - ang,
        horizontal_angles=np.array([0]),
        data=data[:, None],
    )


def to_ascii(filename, apd, /, **kwargs):
    ang = np.arange(181)
    data = np.interp(ang, apd.vertical_angles, apd.vertical_profile(), left=0, right=0)
    np.savetxt(filename, np.stack((data, 180 - ang), axis=1), **kwargs)


def cycle(apd, /):
    ha = apd.horizontal_angles
    data = apd.data

    if len(ha) == 1:
        ha = np.array([0, 90])
        data = np.repeat(data, 2, axis=1)

    rad = np.deg2rad(ha)
    sin = np.sin(rad)
    cos = np.cos(rad)

    if np.abs(np.min(cos)) < 1e-5 or np.max(cos) < 1e-5:
        ha = np.concatenate([ha[:-1], 180 - ha[::-1]])
        data = np.concatenate([data[:, :-1], data[:, ::-1]], axis=1)
    if np.abs(np.min(sin)) < 1e-5 or np.max(sin) < 1e-5:
        ha = np.concatenate([ha[:-1], 360 - ha[::-1]])
        data = np.concatenate([data[:, :-1], data[:, ::-1]], axis=1)

    ha %= 360
    idx = np.argsort(ha)
    data = data[:, idx]
    data = np.concatenate([data, data[:, 0:1]], axis=1)
    ha = ha[idx]
    ha = np.concatenate([ha, ha[0:1] + 360])
    ha, idx = np.unique(ha, return_index=True)
    data = data[:, idx]

    return AngularPowerDistribution(
        lumens=apd.lumens,
        vertical_angles=apd.vertical_angles,
        horizontal_angles=ha,
        data=data,
    )


def vertical_profile(apd, /, *, integrated=False):
    profile = (
        np.average(apd.data, axis=1, weights=np.diff(mids(apd.horizontal_angles)))
        if len(apd.horizontal_angles) > 1
        else apd.data[:, 0].copy()
    )
    if integrated:
        profile *= 2 * np.pi * np.diff(-np.cos(np.deg2rad(mids(apd.vertical_angles))))
    return profile


def horizontal_profile(apd, /):
    profile = apd.data.mean(0)
    return profile / profile.max()


def normalize(apd):
    data = apd.data / np.sum(apd.vertical_profile(integrated=True))

    return AngularPowerDistribution(
        type=apd.type,
        vertical_angles=apd.vertical_angles,
        horizontal_angles=apd.horizontal_angles,
        data=data,
    )


def interpolate(apd, /, *, step=1, method="linear"):
    N = round(90 / step)
    apd = apd.cycle()
    ha, va = np.meshgrid(apd.horizontal_angles, apd.vertical_angles)
    H = np.linspace(0, 360, 4 * N + 1)
    V = np.linspace(0, 180, 2 * N + 1)

    interp = scipy.interpolate.griddata(
        (ha.flatten(), va.flatten()),
        apd.data.flatten(),
        tuple(np.meshgrid(H, V)),
        method=method,
        fill_value=0,
    )

    return AngularPowerDistribution(
        lumens=apd.lumens, vertical_angles=V, horizontal_angles=H, data=interp
    )


def plotV(apd, /, polar=False, **kwargs):
    if polar:
        mpl.pyplot.polar(
            np.deg2rad(apd.vertical_angles), apd.vertical_profile(), **kwargs
        )
    else:
        mpl.pyplot.plot(apd.vertical_angles, apd.vertical_profile(), **kwargs)


def plotH(apd, /, polar=False, **kwargs):
    if polar:
        mpl.pyplot.polar(
            np.deg2rad(apd.horizontal_angles), apd.horizontal_profile(), **kwargs
        )
    else:
        mpl.pyplot.plot(apd.horizontal_angles, apd.horizontal_profile(), **kwargs)


def plot2d(
    apd,
    /,
    ax=None,
    *,
    wrap=False,
    interpolation="nearest",
    cmap=mpl.rc_params()["image.cmap"],
    vmin=0,
    vmax=None,
):
    ha = apd.horizontal_angles
    va = apd.vertical_angles
    data = apd.data

    if wrap and ha[-1] == 360:
        ha = np.roll(ha[:-1], len(ha) // 2)
        ha = (ha + 180) % 360 - 180
        ha = np.concatenate([ha, ha[0:1] + 360])
        data = np.roll(data[:, :-1], len(ha) // 2, axis=1)
        data = np.concatenate([data, data[:, 0:1]], axis=1)

    if ax is None:
        ax = mpl.pyplot.gca()

    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    im = mpl.image.NonUniformImage(
        ax,
        cmap=cmap,
        norm=norm,
        interpolation=interpolation,
        extent=[ha[0], ha[-1], va[0], va[-1]],
    )

    im.set_data(ha, va, data)
    ax.add_image(im)
    ax.set_xlim(ha[0], ha[-1])
    ax.set_ylim(va[0], va[-1])

    return im


def plot3d(apd, /, *, wireframe=False, **kwargs):
    ha, va = np.meshgrid(
        np.deg2rad(90 - apd.horizontal_angles),
        np.deg2rad(180 - apd.vertical_angles),
        copy=False,
        sparse=True,
    )
    r = apd.data
    x = r * np.sin(va) * np.cos(ha)
    y = r * np.sin(va) * np.sin(ha)
    z = r * np.cos(va)

    m = max(np.max(np.abs(x)), np.max(np.abs(y)), np.max(np.abs(z)))

    fig = mpl.pyplot.figure()
    ax = fig.add_subplot(projection="3d")
    (ax.plot_wireframe if wireframe else ax.plot_surface)(
        x, y, z, rstride=1, cstride=1, **kwargs
    )
    ax.set_xlim(-m, m)
    ax.set_ylim(-m, m)
    ax.set_zlim(-m, m)
