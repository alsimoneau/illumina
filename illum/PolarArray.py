#!/usr/bin/env python3

from dataclasses import dataclass

import numpy as np


@dataclass
class PolarArray:
    maxRadius: float
    crs: object
    transform: object
    center: tuple
    shape: tuple
    data: np.array
