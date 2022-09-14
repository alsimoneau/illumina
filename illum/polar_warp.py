#!/usr/bin/env python3

import cv2
import illum.compute
import numpy as np


def warp_polar(image, outshape, center, rmax, *, log=True, inverse=False):
    return cv2.warpPolar(
        src=image,
        dsize=outshape[::-1],  # rho, phi or x y
        center=center[::-1],
        maxRadius=rmax,
        flags=(cv2.WARP_POLAR_LOG if log else cv2.WARP_POLAR_LINEAR)
        | (cv2.WARP_INVERSE_MAP if inverse else 0)
        | cv2.WARP_FILL_OUTLIERS,
    )


def blur_polar(image, outshape, center, rmax, log=True):
    indices = np.arange(1, 1 + np.product(outshape)).reshape(outshape)
    warped_idx = (
        warp_polar(
            indices, image.shape[:2], center, rmax, log=log, inverse=True
        )
        + 1
    )

    data = image.astype("float64")
    if image.ndim > 2:
        avg_image = np.stack(
            [
                illum.compute.average_index(
                    1 + np.product(outshape), warped_idx.T, data[:, :, idx].T
                ).T
                for idx in np.ndindex(image.shape[2:])
            ],
            axis=-1,
        ).reshape(image.shape)
    else:
        avg_image = illum.compute.average_index(
            1 + np.product(outshape), warped_idx.T, data.T
        ).T

    return avg_image.astype(image.dtype)


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    image = np.random.random(100, 100)
    R, T = 1000, 1000
    rmax = min(image.shape[0], image.shape[1]) / 2
    center = (image.shape[0] - 1) / 2, (image.shape[1] - 1) / 2

    image_cv = warp_polar(image, (T, R), center, rmax)
    inage_cv_inv = warp_polar(
        image_cv, image.shape[:2], center, rmax, inverse=True
    )
    image_blurred = blur_polar(image, (T, R), center, rmax)
    image_warped = warp_polar(image_blurred, (T, R), center, rmax)
    image_warped_inv = warp_polar(
        image_warped, image.shape[:2], center, rmax, inverse=True
    )

    fig, axes = plt.subplots(2, 3)
    axes[0, 0].set_title("Original")
    axes[0, 0].imshow(image)
    axes[0, 1].set_title("Polar")
    axes[0, 1].imshow(image_cv)
    axes[0, 2].set_title("Inverse")
    axes[0, 2].imshow(inage_cv_inv)
    axes[1, 0].set_title("Blurred")
    axes[1, 0].imshow(image_blurred)
    axes[1, 1].set_title("Polar")
    axes[1, 1].imshow(image_warped)
    axes[1, 2].set_title("Inverse")
    axes[1, 2].imshow(image_warped_inv)
