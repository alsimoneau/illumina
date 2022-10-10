#!/usr/bin/env python3

import numpy as np

from . import compute


def average_index(arr, indices):
    idx = np.asarray(indices, int)
    n = int(np.prod(arr.shape[idx.ndim :]))
    return (
        compute.average_index(arr.reshape((-1, n)).T, idx.flatten(), idx.max())
        .T.reshape(arr.shape)
        .astype(arr.dtype)
    )


def viirs2lum(*args, **kwargs):
    return compute.viirs2lum(*args, **kwargs)
