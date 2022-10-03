#!/usr/bin/env python3

import numpy as np


def cart2pol(y, x):
    t = np.arctan2(x, y)
    r = np.hypot(x, y)
    return t, r


def pol2cart(t, r):
    y = r * np.cos(t)
    x = r * np.sin(t)
    return y, x


def idx2cart(i, j, dy, dx):
    y = (i + 0.5) * dy
    x = (j + 0.5) * dx
    return y, x


def cart2idx(y, x, dy, dx):
    i = y / dy - 0.5
    j = x / dx - 0.5
    return i, j


def idx2pol(i, j, Nt, Nr, rmin, rmax):
    t = i * 2 * np.pi / Nt
    r = rmin * np.power(rmax / rmin, (j + 0.5) / Nr)
    return t, r


def pol2idx(t, r, Nt, Nr, rmin, rmax):
    i = Nt * t / (2 * np.pi)
    j = Nr * np.log(r / rmin) / np.log(rmax / rmin) - 0.5
    return i, j
