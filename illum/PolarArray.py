#!/usr/bin/env python3

from dataclasses import dataclass

import cv2
import numpy as np

import illum.compute


@dataclass
class PolarArray:
    maxRadius: float
    center: tuple
    shape: tuple
    data: np.array

    def to_array(self):
        return to_array(self)

    def coords(self):
        return coords(self)


def from_array(arr, /, outshape, *, center=None, rmax=None):
    center, rmax = _polar_defaults(arr.shape, center, rmax)
    blur = _blur_polar(arr, outshape, center=center, rmax=rmax, log=True)
    data = _warp_polar(blur, outshape, center=center, rmax=rmax, log=True)
    return PolarArray(
        data=data,
        maxRadius=rmax,
        center=center,
        shape=arr.shape,
    )


def to_array(parr):
    return _warp_polar(
        parr.data, parr.shape, parr.center, parr.rmax, log=True, inverse=True
    )


def coords(parr):
    return _map_coordinates(
        parr.data.shape,
        parr.shape,
        center=parr.center,
        rmax=parr.rmax,
        log=True,
    )


# Utility functions


def _polar_defaults(shape, center=None, rmax=None):
    shape2d = np.array(shape[:2])
    if center is None:
        center = (shape2d - 1) / 2
    if rmax is None:
        rmax = np.sqrt(
            np.sum((np.max((center, shape2d - 1 - center), axis=0) + 0.5) ** 2)
        )

    return np.array(center), rmax


def _correct_rmax(rmax, shape):
    return np.power(rmax + 1, shape[1] / (shape[1] - 0.5))


def _map_coordinates(
    inshape, outshape, /, *, center=None, rmax=None, log=True, inverse=False
):
    if inverse:
        center, rmax = _polar_defaults(outshape, center, rmax)
        rmax = _correct_rmax(rmax, inshape)

        y = np.arange(outshape[0])[:, None] - center[0]
        x = np.arange(outshape[1]) - center[1]

        theta = np.arctan2(y, x) % (2 * np.pi) * inshape[0] / (2 * np.pi)
        r = np.sqrt(x * x + y * y)
        if log:
            r = np.log(r + 1) * inshape[1] / np.log(rmax - 1)
        else:
            r = r * inshape[1] / rmax

        return np.array([theta, r])

    else:
        center, rmax = _polar_defaults(inshape, center, rmax)
        rmax = _correct_rmax(rmax, outshape)

        theta = np.linspace(0, 2 * np.pi, outshape[0], endpoint=False)[:, None]
        if log:
            r = np.geomspace(1, rmax, outshape[1], endpoint=False) - 1
        else:
            r = np.linspace(0, rmax, outshape[1], endpoint=False)

        y = -r * np.sin(theta) + center[0]
        x = r * np.cos(theta) + center[1]

        return np.array([y, x])


def _warp_polar(image, outshape, center, rmax, *, log=True, inverse=False):
    rmax = _correct_rmax(rmax, image.shape if inverse else outshape)
    return cv2.warpPolar(
        src=image,
        dsize=outshape[::-1],  # rho, phi or x y
        center=center[::-1],
        maxRadius=rmax,
        flags=(cv2.WARP_POLAR_LOG if log else cv2.WARP_POLAR_LINEAR)
        | (cv2.WARP_INVERSE_MAP if inverse else 0)
        | cv2.WARP_FILL_OUTLIERS,
    )


def _blur_polar(image, outshape, center, rmax, log=True):
    indices = 1 + np.arange(np.prod(outshape)).reshape(outshape)
    warped_idx = _warp_polar(
        indices, image.shape[:2], center, rmax, log=log, inverse=True
    )

    avg_image = illum.compute.average_index(
        image.reshape((-1, np.prod(image.shape[2:], dtype="uint8"))).T,
        warped_idx.flatten(),
        np.prod(outshape) + 1,
    ).T.reshape(image.shape)

    return avg_image.astype(image.dtype)


def _test_plot(image=None, center=None, rmax=None, res=100):
    import matplotlib.pyplot as plt

    if image is None:
        image = np.random.random((1000, 1000))

    center, rmax = _polar_defaults(image.shape, center, rmax)
    shape = (res, res)

    image_cv = from_array(image, shape, center, rmax)
    inage_cv_inv = to_array(image_cv)
    image_blurred = _blur_polar(image, shape, center, rmax, log=True)
    image_warped = from_array(image_blurred, shape, center, rmax)
    image_warped_inv = to_array(image_warped)

    fig, axes = plt.subplots(2, 3)
    axes[0, 0].set_title("Original")
    axes[0, 0].imshow(image)
    axes[0, 0].plot(*center[::-1], "Pk")
    axes[0, 0].plot(*center[::-1], "+w")
    axes[0, 1].set_title("Polar")
    axes[0, 1].imshow(image_cv.data)
    axes[0, 2].set_title("Inverse")
    axes[0, 2].imshow(inage_cv_inv)
    axes[0, 2].plot(*center[::-1], "Pk")
    axes[0, 2].plot(*center[::-1], "+w")
    axes[1, 0].set_title("Blurred")
    axes[1, 0].imshow(image_blurred)
    axes[1, 0].plot(*center[::-1], "Pk")
    axes[1, 0].plot(*center[::-1], "+w")
    axes[1, 1].set_title("Polar")
    axes[1, 1].imshow(image_warped.data)
    axes[1, 2].set_title("Inverse")
    axes[1, 2].imshow(image_warped_inv)
    axes[1, 2].plot(*center[::-1], "Pk")
    axes[1, 2].plot(*center[::-1], "+w")


if __name__ == "__main__":
    _test_plot()
