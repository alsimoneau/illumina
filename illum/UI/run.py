#!/usr/bin/env python3

import illum
from illum.utils.utils import parallelize


@parallelize
def run(input):
    print(input)
    illum.compute.kernel()
