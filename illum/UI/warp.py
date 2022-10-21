#!/usr/bin/env python3

import os

import fiona
import geopandas as gpd
import numpy as np
import rasterio as rio
import rasterio.features
import rasterio.merge
import rasterio.warp
import shapely.geometry

import illum

fiona.drvsupport.supported_drivers["KML"] = "rw"


def merge(name):
    files = os.listdir(name)
    arr, transform = rio.merge.merge(
        [rio.open(os.path.join(name, s)) for s in files], nodata=0
    )
    crs = rio.open(os.path.join(name, files[0])).crs
    return arr[0], transform, crs


def open_multiple(path):
    outs = []

    A = None
    for name in os.listdir(path):
        name = os.path.join(path, name)
        if os.path.isfile(name):
            rst = rio.open(name)
            if A is None:
                A = rst.transform.determinant
                continue
            if A != rst.transform.determinant:
                break
        else:
            break
    else:
        outs.append(merge(name))
        return outs

    for name in os.listdir(path):
        name = os.path.join(path, name)
        if os.path.isfile(name):
            rst = rio.open(name)
            outs.append((rst.read(1), rst.transform, rst.crs))
        else:
            outs.append(merge(name))
    return outs


def reproject(
    source,
    src_transform,
    src_crs,
    bounds,
    src_nodata=None,
    dst_crs=None,
    dst_nodata=None,
    resampling="average",
):
    bounds = rio.warp.transform_bounds(dst_crs, src_crs, *bounds)
    window = (
        rio.windows.from_bounds(*bounds, src_transform)
        .round_offsets()
        .round_lengths()
    )
    transform, w, h = rio.warp.calculate_default_transform(
        src_crs, dst_crs, window.width, window.height, *bounds
    )
    if not rio.windows.intersect(window, rio.windows.get_data_window(source)):
        return None, None, None

    arr = np.zeros((h, w))
    rio.warp.reproject(
        source[window.toslices()],
        arr,
        src_transform=rio.windows.transform(window, src_transform),
        src_crs=src_crs,
        src_nodata=src_nodata,
        dst_transform=transform,
        dst_crs=dst_crs,
        dst_nodata=dst_nodata,
        resampling=rio.enums.Resampling[resampling],
    )

    if src_nodata is None:
        mask = np.ones(source.shape, dtype=bool)
    elif np.isnan(src_nodata):
        mask = ~np.isnan(source)
    else:
        mask = source != src_nodata

    out_mask = np.zeros((h, w), dtype="uint8")
    rio.warp.reproject(
        mask.astype("uint8"),
        out_mask,
        src_transform=src_transform,
        src_crs=src_crs,
        dst_transform=transform,
        dst_crs=dst_crs,
        resampling=rio.enums.Resampling.nearest,
    )

    return arr, transform, out_mask


def burn(polygon, polarArray, all_touched=False, N=1000, sf=2):
    radii = polarArray.radii()
    area = polarArray.area()
    center = polarArray.center

    outs = []
    i = 0
    while True:
        res = np.sqrt(area[i] / sf)
        transform = rio.transform.from_origin(
            center[1] - N * res, center[0] + N * res, res, res
        )

        arr = rio.features.rasterize(
            polygon,
            (2 * N, 2 * N),
            transform=transform,
            all_touched=all_touched,
        )
        mask = np.ones((2 * N, 2 * N), dtype="bool")
        outs.append((arr, transform, mask))

        try:
            i = np.where(radii > N * res)[0][0]
        except IndexError:
            break

    return outs


def open_polygon(filename, shape, transform, crs):
    mask = gpd.GeoSeries(
        shapely.geometry.Polygon.from_bounds(
            *rio.transform.array_bounds(*shape, transform)
        ),
        crs=crs,
    )
    poly = gpd.read_file(filename, mask=mask).to_crs(crs)
    return poly, mask


def process_water(poly, mask, crs):
    if len(poly) == 0:
        return mask
    elif len(poly) > 1:
        poly = gpd.GeoSeries(poly.unary_union, crs=crs)
    return mask.difference(poly)


def polarize(path, parr, field=None):
    if os.path.isdir(path):
        # multiple rasters
        srcs = [
            reproject(*src, bounds=parr.bounds(), dst_crs=parr.crs)
            for src in open_multiple(path)
        ]
        return illum.PA.union(*zip(*srcs), sort=True, out=parr.copy())

    try:
        rst = rio.open(path)
    except rio.RasterioIOError:
        # vector
        polygon, mask = open_polygon(
            path, parr.xyshape, parr.transform, parr.crs
        )
        if field is None:
            polygon = process_water(polygon, mask, parr.crs)
        else:
            polygon = list(zip(polygon.geometry, polygon[field]))

        srcs = burn(polygon, parr, all_touched=field is None)
        return illum.PA.union(*zip(*srcs), out=parr.copy())

    # raster
    arr, transform, center = reproject(
        rst.read(1),
        rst.transform,
        rst.crs,
        bounds=parr.bounds(),
        dst_crs=parr.crs,
    )
    return illum.PA.from_array(
        arr,
        parr.shape,
        center=parr.center,
        rmin=parr.minRadius,
        rmax=parr.maxRadius,
        crs=parr.crs,
        transform=transform,
    )


def warp(infiles, output=None, domain="domain.parr", field=None):
    polarArray = (
        illum.PA.load(domain)
        if domain.rpartition(".")[2] == "parr"
        else illum.domain(domain)
    )
    parr = polarize(infiles, polarArray, field)
    if output is not None:
        parr.save(output)
    return parr
