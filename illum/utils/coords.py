#!/usr/bin/env python3

import numpy as _np
from scipy.ndimage import map_coordinates as _map_coordinates

import illum


def cart2pol(coords, center=(0, 0)):
    y = coords[0] - center[0]
    x = coords[1] - center[1]
    t = _np.arctan2(x, y) % (2 * _np.pi)
    r = _np.hypot(x, y)
    return t, r


def pol2cart(coords, center=(0, 0)):
    y = coords[1] * _np.cos(coords[0])
    x = coords[1] * _np.sin(coords[0])
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
    t = idx[0] * 2 * _np.pi / shape[0]
    r = lims[0] * _np.power(lims[1] / lims[0], (idx[1] + 0.5) / shape[1])
    return t, r


def pol2idx(coords, shape, lims):
    i = shape[0] * coords[0] / (2 * _np.pi)
    dr = _np.log(coords[1] / lims[0]) / _np.log(lims[1] / lims[0])
    j = shape[1] * dr - 0.5
    return i, j


def wrap(arr, shape, center, origin, scale, lims, res=None, order=1):
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
    idx = _np.indices(shape, sparse=True)
    coords = pol2cart(idx2pol(idx, res, lims), center)
    idx = cart2idx(coords, origin, scale)
    if order == 0:
        idx = _np.round(idx)

    return _map_coordinates(arr, idx, order=order)


def unwrap(arr, shape, center, origin, scale, lims, res=None, order=1):
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
    idx = _np.indices(shape, sparse=True)
    coords = cart2pol(idx2cart(idx, origin, scale), center)
    idx = pol2idx(coords, res, lims)
    if order == 0:
        idx = _np.round(idx)

    arr = _np.vstack((arr, arr[0:1]))  # Polar cycle
    return _map_coordinates(arr, idx, order=order)


def polar_blur(arr, shape, center, origin, scale, lims, res=None):
    idx = 1 + _np.arange(_np.prod(shape)).reshape(shape)
    w_idx = unwrap(idx, arr.shape, center, origin, scale, lims, res, order=0)

    avg = illum.compute.average_index(arr, w_idx)
    avg[w_idx == 0] = 0
    return avg
