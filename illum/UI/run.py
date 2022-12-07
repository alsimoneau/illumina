#!/usr/bin/env python3

from illum.compute import kernel
from illum.utils import parallelize


@parallelize
def run(input):
    print(input)
