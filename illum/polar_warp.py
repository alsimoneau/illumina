#!/usr/bin/env python3

import cv2
import numpy as np
import scipy.ndimage

import illum.compute


def warp_polar(image, R, T, center, rmax, *, log=True, inverse=False):
    return cv2.warpPolar(
        src=image,
        dsize=(R, T),  # rho, phi
        center=center,
        maxRadius=rmax,
        flags=(cv2.WARP_POLAR_LOG if log else cv2.WARP_POLAR_LINEAR)
        | (cv2.WARP_INVERSE_MAP if inverse else 0)
        | cv2.WARP_FILL_OUTLIERS,
    )


def blur_polar(image, R, T, center, rmax, log=True):
    indices = np.arange(1, 1 + R * T).reshape((R, T))
    warped_idx = (
        warp_cv(indices, *image.shape[:2], center, rmax, log=log, inverse=True)
        + 1
    )

    data = image.astype("float64")
    if image.ndim > 2:
        avg_image = np.stack(
            [
                compute.average_index(
                    R * T + 1, warped_idx.T, data[:, :, idx].T
                )
                for idx in np.ndindex(image.shape[2:])
            ],
            axis=-1,
        ).reshape(image.shape)
    else:
        avg_image = illum.compute.average_index(
            R * T + 1, warped_idx.T, data.T
        )

    return avg_image.astype(image.dtype)


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    image = np.random.random(100, 100)
    R, T = 1000, 1000
    rmax = image.shape[1] / 2
    center = (image.shape[0] - 1) / 2, (image.shape[1] - 1) / 2

    image_cv = warp_polar(image, R, T, center, rmax)
    image_blurred = blur_polar(image, R, T, center, rmax)
    image_warped = warp_polar(image_blurred, R, T, center, rmax)

    fig, axes = plt.subplots(2, 2)
    axes[0, 0].set_title("Original")
    axes[0, 0].imshow(image)
    axes[0, 1].set_title("OpenCV")
    axes[0, 1].imshow(image_cv)
    axes[1, 0].set_title("Blurred")
    axes[1, 0].imshow(image_blurred)
    axes[1, 1].set_title("Warped")
    axes[1, 1].imshow(image_warped)
