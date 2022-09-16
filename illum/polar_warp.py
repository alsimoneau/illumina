#!/usr/bin/env python3

import cv2
import numpy as np

import illum.compute


def warp_polar(image, outshape, center, rmax, *, log=True, inverse=False):
    rmax = correct_rmax(rmax, image.shape if inverse else outshape)
    return cv2.warpPolar(
        src=image,
        dsize=outshape[::-1],  # rho, phi or x y
        center=center[::-1],
        maxRadius=rmax,
        flags=(cv2.WARP_POLAR_LOG if log else cv2.WARP_POLAR_LINEAR)
        | (cv2.WARP_INVERSE_MAP if inverse else 0)
        | cv2.WARP_FILL_OUTLIERS,
    )


def map_coordinates(
    inshape, outshape, /, *, center=None, rmax=None, log=True, inverse=False
):
    if inverse:
        center, rmax = polar_defaults(outshape, center, rmax)
        rmax = correct_rmax(rmax, inshape)

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
        center, rmax = polar_defaults(inshape, center, rmax)
        rmax = correct_rmax(rmax, outshape)

        theta = np.linspace(0, 2 * np.pi, outshape[0], endpoint=False)[:, None]
        if log:
            r = np.geomspace(1, rmax, outshape[1], endpoint=False) - 1
        else:
            r = np.linspace(0, rmax, outshape[1], endpoint=False)

        y = -r * np.sin(theta) + center[0]
        x = r * np.cos(theta) + center[1]

        return np.array([y, x])


def blur_polar(image, outshape, center, rmax, log=True):
    indices = 1 + np.arange(np.product(outshape)).reshape(outshape)
    warped_idx = warp_polar(
        indices, image.shape[:2], center, rmax, log=log, inverse=True
    )

    if image.ndim > 2:
        avg_image = np.stack(
            [
                illum.compute.average_index(
                    np.product(outshape) + 1, warped_idx.T, image[:, :, idx].T
                ).T
                for idx in np.ndindex(image.shape[2:])
            ],
            axis=-1,
        ).reshape(image.shape)
    else:
        avg_image = illum.compute.average_index(
            np.product(outshape) + 1, warped_idx.T, image.T
        ).T

    return avg_image.astype(image.dtype)


def polar_defaults(shape, center=None, rmax=None):
    shape2d = np.array(shape[:2])
    if center is None:
        center = (shape2d - 1) / 2
    if rmax is None:
        rmax = np.sqrt(
            np.sum((np.max((center, shape2d - 1 - center), axis=0) + 0.5) ** 2)
        )

    return np.array(center), rmax


def correct_rmax(rmax, shape):
    return np.power(rmax + 1, shape[1] / (shape[1] - 0.5))


def polar_warp(image, /, shape, *, center=None, rmax=None, log=True):
    center, rmax = polar_defaults(image.shape, center, rmax)
    blurred = blur_polar(image, shape, center, rmax, log=log)
    return warp_polar(blurred, shape, center, rmax, log=log)


def polar_unwarp(image, /, shape, *, center=None, rmax=None, log=True):
    center, rmax = polar_defaults(shape, center, rmax)
    return warp_polar(image, shape, center, rmax, log=log, inverse=True)


def plot_test(image=None, center=None, rmax=None, res=100, log=True):
    import matplotlib.pyplot as plt

    if image is None:
        image = np.random.random((1000, 1000))
    center, rmax = polar_defaults(image.shape, center, rmax)
    shape = (res, res)

    image_cv = warp_polar(image, shape, center, rmax, log=log)
    inage_cv_inv = warp_polar(
        image_cv, image.shape[:2], center, rmax, log=log, inverse=True
    )
    image_blurred = blur_polar(image, shape, center, rmax, log=log)
    image_warped = warp_polar(image_blurred, shape, center, rmax, log=log)
    image_warped_inv = warp_polar(
        image_warped, image.shape[:2], center, rmax, log=log, inverse=True
    )

    fig, axes = plt.subplots(2, 3)
    axes[0, 0].set_title("Original")
    axes[0, 0].imshow(image)
    axes[0, 0].plot(*center[::-1], "Pk")
    axes[0, 0].plot(*center[::-1], "+w")
    axes[0, 1].set_title("Polar")
    axes[0, 1].imshow(image_cv)
    axes[0, 2].set_title("Inverse")
    axes[0, 2].imshow(inage_cv_inv)
    axes[0, 2].plot(*center[::-1], "Pk")
    axes[0, 2].plot(*center[::-1], "+w")
    axes[1, 0].set_title("Blurred")
    axes[1, 0].imshow(image_blurred)
    axes[1, 0].plot(*center[::-1], "Pk")
    axes[1, 0].plot(*center[::-1], "+w")
    axes[1, 1].set_title("Polar")
    axes[1, 1].imshow(image_warped)
    axes[1, 2].set_title("Inverse")
    axes[1, 2].imshow(image_warped_inv)
    axes[1, 2].plot(*center[::-1], "Pk")
    axes[1, 2].plot(*center[::-1], "+w")


if __name__ == "__main__":
    plot_test()
