# General convenience functions
# Author: Alexandre Simoneau


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
