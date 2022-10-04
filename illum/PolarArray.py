#!/usr/bin/env python3

from copy import copy
from dataclasses import dataclass

import h5py
import numpy as np
import rasterio as rio

import illum


@dataclass
class PolarArray:
    data: np.ndarray
    xyshape: tuple
    center: tuple
    minRadius: float
    maxRadius: float
    crs: object = "None"
    transform: object = None

    def __post_init__(self):
        if self.transform is None:
            self.transform = rio.transform.IDENTITY
        else:
            try:
                self.transform = rio.transform.Affine(*self.transform)
            except TypeError:
                pass

    def __repr__(self):
        return (
            f"PolarArray<scale={self._scale},"
            f"[{self.shape[0]},{self.minRadius}:{self.maxRadius}]>"
        )

    @property
    def shape(self):
        return self.data.shape

    @property
    def _origin(self):
        return self.transform.f, self.transform.c

    @property
    def _scale(self):
        return self.transform.e, self.transform.a

    @property
    def _lims(self):
        return self.minRadius, self.maxRadius

    @property
    def center_coord(self):
        return illum.utils.cart2idx(self.center, self._origin, self._scale)

    def radii(self):
        return np.geomspace(
            self.minRadius, self.maxRadius, 2 * self.shape[1] + 1
        )[1::2]

    def area(self):
        r = np.geomspace(self.minRadius, self.maxRadius, self.shape[1] + 1)
        return (
            (np.pi / self.shape[0])
            * np.diff(r ** 2)
            * np.abs(np.prod(self._scale))
        )

    def copy(self):
        return copy(self)

    def to_array(self):
        return to_array(self)

    def save(self, filename, *args, **kwargs):
        return save(filename, self, *args, **kwargs)


def load(filename):
    with h5py.File(filename, "r") as f:
        parr = PolarArray(data=f["data"][:], **f.attrs)
    return parr


def save(filename, parr):
    if filename[:-5] != ".parr":
        filename += ".parr"

    with h5py.File(filename, "w") as f:
        attrs = parr.__dict__
        f.create_dataset("data", data=attrs["data"])
        if attrs["transform"] != "None":
            attrs["transform"] = tuple(attrs["transform"])[:-3]
        attrs = {k: v for k, v in attrs.items() if k[0] != "_" and k != "data"}
        f.attrs.update(attrs)


def from_array(
    arr,
    /,
    outshape,
    *,
    center=None,
    rmin=0,
    rmax=None,
    crs="None",
    transform=None,
):
    if type(outshape) == PolarArray:
        rmin = outshape.minRadius
        rmax = outshape.maxRadius
        center = outshape.center
        transform = outshape.transform
        crs = outshape.crs

    if transform is None:
        transform = rio.transform.IDENTITY

    origin = (transform.f, transform.c)
    scale = (transform.e, transform.a)
    lims = (rmin, rmax)

    blur = illum.utils.polar_blur(arr, outshape, center, origin, scale, lims)
    data = illum.utils.wrap(blur, outshape, center, origin, scale, lims)

    return PolarArray(
        data=data,
        minRadius=rmin,
        maxRadius=rmax,
        center=center,
        xyshape=arr.shape,
        crs=crs,
        transform=transform,
    )


def to_array(parr):
    return illum.utils.unwrap(
        parr.data,
        parr.xyshape,
        parr.center,
        parr._origin,
        parr._scale,
        parr._lims,
    )


def union(
    arrays,
    /,
    transforms=None,
    centers=None,
    masks=None,
    *,
    sort=False,
    out=None,
    shape=None,
    maxRadius=None,
    crs=None,
):
    if not (transforms is None and centers is None) and (
        len(arrays) != len(transforms) != len(centers)
    ):
        raise ValueError(
            "'arrays', 'transforms' and 'centers' must have the same length."
        )

    if out is None:
        if not (shape is None or maxRadius is None or crs is None):
            raise ValueError(
                "'shape', 'maxRadius' and 'crs' must be given when 'out' is"
                " not provided."
            )
        else:
            out = from_array([[0]], shape, rmax=maxRadius, crs=crs)
    else:
        if not isinstance(out, PolarArray):
            raise TypeError("'out' must be a PolarArray.")

    if masks is None:
        masks = [None] * len(arrays)

    if sort:
        areas = [np.abs(t.determinant) for t in transforms]
        idx = np.argsort(areas)

        arrays = arrays[idx]
        masks = masks[idx]
        transforms = transforms[idx]

    for arr, mask, center, transform in zip(
        arrays[::-1], masks[::-1], centers[::-1], transforms[::-1]
    ):
        if mask is None:
            mask = np.ones(arr.shape)

        domain = from_array(
            mask,
            out.shape,
            center=center,
            rmax=out.maxRadius * out._scale,
            crs=out.crs,
            transform=transform,
        )
        print(domain.radii() / domain._scale)
        domain = domain.data
        data = from_array(
            arr,
            out.shape,
            center=center,
            rmax=out.maxRadius * out._scale,
            crs=out.crs,
            transform=transform,
        ).data
        pmask = domain > 0.5
        out.data[pmask] = data[pmask]

    return out


def _test_polar():
    new = from_array(
        np.random.random((100, 100)),
        (100, 100),
        center=(50, 50),
        rmin=1,
        rmax=50,
    )
    new.save("toto")
    old = load("toto.parr")
    np.allclose(old.data, new.data)


if __name__ == "__main__":
    _test_polar()
