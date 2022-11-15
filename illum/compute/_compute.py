#!/usr/bin/env python3

import numpy as np

from .compute import viirs2lum


def average_index(arr, indices):
    from .compute import average_index

    idx = np.asarray(indices, int)
    n = int(np.prod(arr.shape[idx.ndim :]))
    return (
        average_index(arr.reshape((-1, n)).T, idx.flatten(), idx.max())
        .T.reshape(arr.shape)
        .astype(arr.dtype)
    )
