#!/usr/bin/env python3
# ******************************************************************************
#        Determine street orientation from list of coordinates
#
# Authors : Julien-Pierre Houle & Alexandre Simoneau
# Date :    May 2021
# ******************************************************************************

from math import ceil as _ceil

import numpy as _np
import osmnx as _ox
import pyproj as _pyproj
import rasterio as _rio
import rasterio.features
import shapely.geometry as _geometry
import shapely.ops as _ops
from scipy.ndimage import distance_transform_edt as _distance_transform_edt


def transform(s_crs=4326, t_crs=4326):
    s_crs = _pyproj.CRS.from_epsg(s_crs)
    t_crs = _pyproj.CRS.from_epsg(t_crs)
    return _pyproj.Transformer.from_crs(s_crs, t_crs, always_xy=True).transform


# https://gdal.org/tutorials/geotransforms_tut.html
def geotransform(X, Y, GT):
    print("WARNING: Depreciated GDAL geotransform")
    X_geo = GT[0] + (X + 0.5) * GT[1] + (Y + 0.5) * GT[2]
    Y_geo = GT[3] + (X + 0.5) * GT[4] + (Y + 0.5) * GT[5]
    return X_geo, Y_geo


def estimate_utm_epsg(lon, lat):
    code = _pyproj.database.query_utm_crs_info(
        datum_name="WGS84",
        area_of_interest=_pyproj.aoi.AreaOfInterest(
            _np.min(lon), _np.min(lat), _np.max(lon), _np.max(lat)
        ),
    )[0].code
    return int(code)


def fishnet(bounds, res, xsize=None, ysize=None, cols=None, rows=None):
    if xsize is not None and cols is not None:
        print("WARNING: 'xsize' and 'cols' provided. Ignoring 'cols'.")
    if ysize is not None and rows is not None:
        print("WARNING: 'ysize' and 'rows' provided. Ignoring 'rows'.")

    xmin, ymin, xmax, ymax = bounds.bounds
    if cols is None:
        cols = round((xmax - xmin) / xsize)
    xsize = _ceil((xmax - xmin) / cols / res) * res
    if rows is None:
        rows = round((ymax - ymin) / ysize)
    ysize = _ceil((ymax - ymin) / rows / res) * res

    geoms = {}
    for row in range(rows):
        for col in range(cols):
            box = _geometry.box(
                xmin + xsize * col,
                ymax - ysize * (row + 1),
                xmin + xsize * (col + 1),
                ymax - ysize * row,
            )
            geo = bounds.intersection(box)
            if not geo.is_empty:
                geoms[row, col] = [box, geo]

    return geoms, (int(ysize // res), int(xsize // res))


# https://gist.github.com/calebrob6/5039fd409921606aa5843f8fec038c03
def download_roads(bounds, crs):
    bounds_lonlat = _ops.transform(transform(s_crs=crs), bounds)
    Graph = _ox.graph_from_polygon(
        bounds_lonlat,
        network_type="drive",
        simplify=False,
        retain_all=True,
        truncate_by_edge=True,
    )
    return _ox.graph_to_gdfs(Graph, nodes=False)


# https://gis.stackexchange.com/a/151861/92556
def rasterize_roads(roads, bounds, epsg, xres, yres):
    crs = _pyproj.CRS.from_epsg(epsg)
    roads = roads.to_crs(crs)
    xmin, ymin, xmax, ymax = bounds.bounds
    width = int((xmax - xmin) / xres)
    height = int((ymax - ymin) / yres)
    transform = _rio.transform.from_bounds(
        xmin, ymin, xmax, ymax, width, height
    )
    burned = _rio.features.rasterize(
        roads.geometry,
        out_shape=(height, width),
        fill=1,
        default_value=0,
        all_touched=True,
        transform=transform,
    )
    return burned


def dist_ang(roads, xres, yres):
    dist, i = _distance_transform_edt(
        roads, return_indices=True, sampling=(abs(yres), xres)
    )
    Y, X = _np.indices(roads.shape, sparse=True)
    return dist, _np.arctan2((i[1] - X) * xres, (i[0] - Y) * yres)


def compute_ground_type(dist, bounds):
    ground_type = _np.ones(dist.shape) * len(bounds)
    for i, lim in reversed(list(enumerate(sorted(bounds)))):
        ground_type[dist < lim] = i


def nearest_road(lat, lon, radius=300):
    epsg = estimate_utm_epsg(lon, lat)
    x, y = transform(t_crs=epsg)(lon, lat)
    try:
        graph = _ox.graph_from_point(
            (lat, lon),
            dist=radius,
            network_type="drive",
            simplify=False,
            retain_all=True,
            truncate_by_edge=True,
            clean_periphery=True,
        )
    except (ValueError, _ox._errors.EmptyOverpassResponse):
        return _np.nan, 0

    graph = _ox.get_undirected(graph)
    graph = _ox.add_edge_bearings(graph)
    graph = _ox.projection.project_graph(graph, to_crs=epsg)
    edges = _ox.graph_to_gdfs(graph, nodes=False)
    nearest_idx, distance = _ox.distance.nearest_edges(
        graph, x, y, return_dist=True
    )
    nearest = edges.loc[nearest_idx]
    bearing = nearest.bearing
    lon_c, lat_c = transform(s_crs=epsg)(*nearest.geometry.coords[0])
    brng = _ox.bearing.calculate_bearing(lat_c, lon_c, lat, lon)
    if (brng - bearing) % 360 > 180:
        bearing += 180
    bearing -= 90

    return distance, bearing
