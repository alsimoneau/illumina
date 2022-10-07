#!/usr/bin/env python3

import numpy as np
from scipy.ndimage import map_coordinates

import illum


def cart2pol(coords, center=(0, 0)):
    y = coords[0] - center[0]
    x = coords[1] - center[1]
    t = np.arctan2(x, y) % (2 * np.pi)
    r = np.hypot(x, y)
    return t, r


def pol2cart(coords, center=(0, 0)):
    y = coords[1] * np.cos(coords[0])
    x = coords[1] * np.sin(coords[0])
    return y + center[0], x + center[1]


def idx2cart(idx, origin=(0, 0), scale=(-1, 1)):
    y = origin[0] + (idx[0] + 0.5) * scale[0]
    x = origin[1] + (idx[1] + 0.5) * scale[1]
    return y, x


def cart2idx(coords, origin=(0, 0), scale=(-1, 1)):
    i = (coords[0] - origin[0]) / scale[0] - 0.5
    j = (coords[1] - origin[1]) / scale[1] - 0.5
    return i, j


def idx2pol(idx, shape, lims):
    t = idx[0] * 2 * np.pi / shape[0]
    r = lims[0] * np.power(lims[1] / lims[0], (idx[1] + 0.5) / shape[1])
    return t, r


def pol2idx(coords, shape, lims):
    i = shape[0] * coords[0] / (2 * np.pi)
    dr = np.log(coords[1] / lims[0]) / np.log(lims[1] / lims[0])
    j = shape[1] * dr - 0.5
    return i, j


def wrap(arr, shape, center, origin, scale, lims, res=None):
    """
    Wraps data in a polar array.

    Parameters:
        arr : data array
        shape : output shape
        center : polar center point in real coordinates
        origin : cartesian origin
        scale : cartesian scaling
        lims : polar radial limits

    Returns:
        warped array
    """
    if res is None:
        res = shape
    idx = np.indices(shape, sparse=True)
    coords = pol2cart(idx2pol(idx, res, lims), center)
    idx = np.round(cart2idx(coords, origin, scale))

    return map_coordinates(arr, idx, order=0)


def unwrap(arr, shape, center, origin, scale, lims, res=None):
    """
    Unwraps a polar array.

    Parameters:
        arr : data array
        shape : output shape
        center : polar center point in real coordinates
        origin : cartesian origin
        scale : cartesian scaling
        lims : polar radial limits

    Returns:
        warped array
    """
    if res is None:
        res = arr.shape
    idx = np.indices(shape, sparse=True)
    coords = cart2pol(idx2cart(idx, origin, scale), center)
    idx = np.round(pol2idx(coords, res, lims))

    arr = np.vstack((arr, arr[0:1]))  # Polar cycle
    return map_coordinates(arr, idx, order=0)


def polar_blur(arr, shape, center, origin, scale, lims, res=None):
    idx = 1 + np.arange(np.prod(shape)).reshape(shape)
    w_idx = unwrap(idx, arr.shape, center, origin, scale, lims, res)

    avg = illum.compute.average_index(
        arr.reshape((-1, np.prod(arr.shape[2:], dtype="uint32"))).T,
        w_idx.flatten(),
        np.prod(shape) + 1,
    ).T.reshape(arr.shape)

    return avg.astype(arr.dtype)
