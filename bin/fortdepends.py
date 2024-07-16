#!/usr/bin/env python3

import os
import re
import sys
from collections import defaultdict

fnames = sys.argv[1:]

define = re.compile(r"MODULE ([A-Z_0-9]+)")
using = re.compile(r"USE ([A-Z_0-9]+)")

def_names = dict()
use_names = defaultdict(set)

for fname in fnames:
    with open(fname) as f:
        content = f.read().upper()

    for match in re.findall(define, content):
        def_names[match] = fname

    for match in re.findall(using, content):
        use_names[match].add(fname)

use_names = {key: val for key, val in use_names.items() if key in def_names}

depends = {fname: set() for fname in fnames}
for name, fnames in use_names.items():
    for fname in fnames:
        depends[fname].add(name)

depends = {
    fname: {def_names[name] for name in names if def_names[name] != fname}
    for fname, names in depends.items()
}


def libname(fname):
    return f"lib/{os.path.splitext(os.path.basename(fname))[0]}.o"


dep_str = ""
for fname, dep in depends.items():
    dep_str += f"{libname(fname)} :  \\\n"
    for name in dep:
        dep_str += f"\t{libname(name)} \\\n"
    dep_str = dep_str[:-2] + "\n\n"

with open("lib/depends.mk", "w") as f:
    f.write(dep_str)
