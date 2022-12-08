#!/usr/bin/env python3

import numpy as _np


def average_index(arr, indices):
    """Averages an array over like indices.

    Parameters:
        arr: (n+m)D array to average
        indices: nD indices map

    The first n dimensions must be the same.
    """
    from .compute import average_index

    idx = _np.asarray(indices, int)
    n = int(_np.prod(arr.shape[idx.ndim :]))
    return (
        average_index(arr.reshape((-1, n)).T, idx.flatten())
        .T.reshape(arr.shape)
        .astype(arr.dtype)
    )


def viirs2lum(viirs_dat, zones, zonData, inputs):
    """Converts VIIRS data to luminance.

    Parameters:
        viirs_dat: VIIRS data [PolarArray]
        zones: inventory zones [PolarArray]
        zonData: inventory, as returned by 'illum.utils.inputs.parse_inventory'
        inputs: obtained from 'illum.utils.inputs.prep_inputs'
    """
    from .compute import viirs2lum

    spct_keys = list(inputs["spcts"].keys())
    lop_keys = list(inputs["lops"].keys())

    zones_inventory = _np.array(
        [
            [i, lamp[0], spct_keys.index(lamp[1]), lop_keys.index(lamp[2])]
            for i, zon in enumerate(zonData)
            for lamp in zon
        ]
    )
    zones_inventory[:, [0, 2, 3]] += 1  # Fortran indexing

    sources = {lamp[2] for zd in zonData for lamp in zd}
    sources_idx = [list(sources).index(k) if k in sources else -1 for k in lop_keys]
    sources_idx = _np.array(sources_idx) + 1

    return viirs2lum(
        nzones=len(zonData),
        nsources=len(sources),
        viirs=viirs_dat.data,
        zones=zones.data.astype("uint32"),
        angles=inputs["lops"][lop_keys[0]].vertical_angles,
        wav=inputs["spcts"][spct_keys[0]].wavelengths,
        bands=inputs["bins"],
        sens=inputs["sens"].data,
        lops=_np.array([inputs["lops"][k].vertical_profile() for k in lop_keys]),
        spcts=_np.array([inputs["spcts"][k].data for k in spct_keys]),
        sources=sources_idx,
        ivtr=zones_inventory,
        pixsize=zones.area(),
        reflect=inputs["refls"],
    )


def kernel():
    from .compute import kernel
