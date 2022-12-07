# General convenience functions
# Author: Alexandre Simoneau

import math
import os
from glob import glob

import joblib


def strip_comments(item, token="#"):
    """Generator. Strips comments and whitespace from input lines.

    This generator strips comments, leading/trailing whitespace, and
    blank lines from its input.

    Arguments:
        item (obj):  Object to strip comments from.
        token (str, optional):  Comment delimiter.  Defaults to ``#``.

    Yields:
        str:  Next non-blank line from ``item`` with comments and
            leading/trailing whitespace removed.

    Credits: Doug R., StackOverflow
    """

    for line in item:
        s = line.split(token, 1)[0].strip()
        if s != "":
            yield s


def eng_format(x, unit=""):
    # Credit: 200_success on StackOverflow
    # https://codereview.stackexchange.com/a/50971
    #
    # U+03BC is Greek lowercase mu
    UNITS = (
        [" ", " k", " M", " G"]
        + ([None] * 10)
        + [" f", " p", " n", " \u03bc", " m"]
    )

    power_of_1000 = int(math.floor(math.log10(x) // 3))
    exponent = 3 * power_of_1000
    prefix = UNITS[power_of_1000]
    if prefix is None:
        prefix = "*10^%d " % exponent

    significand = x * 10 ** (-exponent)
    return "%.2f%s%s" % (significand, prefix, unit)


def chunker(seq, size):
    return (seq[pos : pos + size] for pos in range(0, len(seq), size))


def add_arrays(a, b):
    if len(a) < len(b):
        a, b = b, a
    c = a.copy()
    c[: len(b)] += b
    return c


def glob_types(p, types):
    return [s for s in glob(p) if os.path.splitext(s)[1][1:].lower() in types]


def parallelize(func):
    def wrapper(iterable, *args):
        return joblib.Parallel(n_jobs=-1, prefer="threads")(
            joblib.delayed(func)(i, *args) for i in iterable
        )

    return wrapper
