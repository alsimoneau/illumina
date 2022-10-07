#!/usr/bin/env python3

import os
from glob import glob

import geopandas as gpd
import numpy as np
import rasterio as rio
import rasterio.features
import rasterio.merge
import rasterio.warp
import shapely.geometry

import illum.PolarArray as PA


def open_multiple(path):
    outs = []
    for name in os.listdir(path):
        name = os.path.join(path, name)
        if os.path.isfile(name):
            rst = rio.open(name)
            outs.append((rst.read(1), rst.transform, rst.crs))
        else:
            files = os.listdir(name)
            arr, transform = rio.merge.merge(
                [rio.open(os.path.join(name, s)) for s in files], nodata=0
            )
            crs = rio.open(os.path.join(name, files[0])).crs
            outs.append((arr[0], transform, crs))
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


def process_water(filename, shape, transform, crs):
    mask = gpd.GeoSeries(
        shapely.geometry.Polygon.from_bounds(
            *rio.transform.array_bounds(*shape, transform)
        ),
        crs=crs,
    )
    poly = gpd.read_file(filename, mask=mask).to_crs(crs)
    if len(poly) == 0:
        return mask
    elif len(poly) > 1:
        poly = gpd.GeoSeries(poly.unary_union, crs=crs)
    return mask.difference(poly)


def polarize(path, parr):
    if os.path.isdir(path):
        # multiple rasters
        srcs = [
            reproject(*src, bounds=parr.bounds(), dst_crs=parr.crs)
            for src in open_multiple(path)
        ]
        return PA.union(*zip(*srcs), sort=True, out=parr.copy())

    try:
        rst = rio.open(path)
    except rio.RasterioIOError:
        # vector
        polygon = process_water(path, parr.xyshape, parr.transform, parr.crs)
        srcs = burn(polygon, parr, all_touched=True)
        return PA.union(*zip(*srcs), out=parr.copy())

    # raster
    arr, transform, center = reproject(
        rst.read(1),
        rst.transform,
        rst.crs,
        bounds=parr.bounds(),
        dst_crs=parr.crs,
    )
    return PA.from_array(
        arr,
        parr.shape,
        center=parr.center,
        rmin=parr.minRadius,
        rmax=parr.maxRadius,
        crs=parr.crs,
        transform=transform,
    )


def correction_filenames(srcfiles):
    try:
        return [
            fname.replace(fname.split("_")[3], "zero_correction").replace(
                "avg_rade9h.tif", "csv"
            )
            for fname in srcfiles
        ]
    except IndexError:
        return [fname.replace("avg_rade9h.tif", "csv") for fname in srcfiles]


def convert_correction_data(srcfiles):
    corr_files = np.unique(correction_filenames(srcfiles))

    data = np.nanmean(
        [np.loadtxt(fname, delimiter=",") for fname in corr_files], 0
    )
    data[np.isnan(data)] = -9999

    with open("VIIRS-DNB/correction.asc", "w") as f:
        f.write(
            "NCOLS 72\n"
            "NROWS 28\n"
            "XLLCORNER -180\n"
            "YLLCORNER -65\n"
            "CELLSIZE 5\n"
            "NODATA_VALUE -9999\n"
        )
        np.savetxt(f, data)

    with open("VIIRS-DNB/correction.prj", "w") as f:
        f.write(
            "GEOGCS["
            '"GCS_WGS_1984",'
            'DATUM["D_WGS_1984",'
            'SPHEROID["WGS_1984",6378137,298.257223563]],'
            'PRIMEM["Greenwich",0],'
            'UNIT["Degree",0.017453292519943295]'
            "]"
        )


def warp(output_name=None, infiles=None):
    if output_name is not None and infiles is None:
        print(
            "ERROR: If an output name is given, files to process must also be"
            " provided."
        )
        raise SystemExit

    polarArray = PA.load("domain.parr")

    if len(infiles):
        polarize(infiles, polarArray).save(output_name)
        return
    else:
        if os.path.isfile("GHSL.zip"):
            print("Found GHSL.zip file, processing.")
            data = [
                warp_files(
                    ["/vsizip/GHSL.zip/GHSL.tif"], params["srs"], extent
                )
                for extent in params["extents"]
            ]
            save(params, data, "obstf")
        else:
            print("WARNING: Could not find GHSL.zip file.")
            print("If you don't intend to use it, you can safely ignore this.")

        files = sorted(glob("SRTM/*.hgt"))
        if not len(files):
            print("ERROR: Could not find SRTM file(s), aborting.")
            raise SystemExit
        print("    ".join(map(str, files)))
        data = [
            warp_files(files, params["srs"], extent)
            for extent in params["extents"]
        ]
        save(params, data, "srtm")

        files = sorted(glob("VIIRS-DNB/*.tif"))
        if not len(files):
            print("WARNING: Did not find VIIRS file(s).")
            print(
                "If you don't intend to use zones inventory, you can safely"
                " ignore this."
            )
        else:
            if not os.path.isfile("hydropolys.zip"):
                print("ERROR: Could not find hydropolys.zip file, aborting.")
                raise SystemExit

            print("    ".join(map(str, files)))

            correction = np.all(
                [
                    os.path.isfile(fname)
                    for fname in correction_filenames(files)
                ]
            )
            if not correction:
                print(
                    "WARNING: Could not find correction files that matched the"
                    " VIIRS files."
                )
                print(
                    "If you indend to use it, please validate that you have"
                    " the right ones."
                )
                print(
                    "Note that only the VCMCFG dataset from VIIRS can be"
                    " corrected."
                )

            data = [
                warp_files(files, params["srs"], extent)
                for extent in params["extents"]
            ]

            if correction:
                convert_correction_data(files)
                corr = [
                    warp_files(
                        ["VIIRS-DNB/correction.asc"], params["srs"], extent
                    )
                    for extent in params["extents"]
                ]
                save(params, corr, "VIIRS_background")
                for i, c in enumerate(corr):
                    data[i] -= c

            save(params, data, "stable_lights")

            prep_shp(
                "hydropolys.zip/hydropolys.shp",
                params["srs"],
                params["extents"][-1],
            )
            data = [
                rasterize("tmp_merge.shp", params["srs"], extent)
                for extent in params["extents"]
            ]
            save(params, data, "water_mask")

        for fname in glob("tmp_*"):
            os.remove(fname)

    print("Done.")
