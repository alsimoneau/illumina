# Functions related to plotting
# Author: Alexandre Simoneau

import matplotlib as _mpl
import matplotlib.colors
import matplotlib.pyplot as _plt
import numpy as _np
from scipy.interpolate import interp2d as _interp2d


def plot_allsky(phi, r, data, n=100, **kwargs):
    """Plot all sky data.

    Parameters:
      phi      : Azimuthal angle (deg). Must be linearly spaced.
      r        : Elevation angle (deg). Must be in [0;90]
      data     : (NxM) or (NxMx3) data to plot with N the number of elevation
                 angles and M the number of azimuthal angles. If 3D the Data
                 will be represented as an RGB plot.
      n        : Supersampling. If 'interp' is not 'None', it will be applied
                 to both axis. (Used to make the plot round)
      interp   : Interpolation method. Can be 'None', 'linear', 'cubic' or
                 'quintic'
      autogain : If True, the data will be stretched to have a max of 1
                 (def. False)
      title    : Figure title
      clabel   : Colorbar label. If defined, will draw a colorbar (2D data
                 only)
      fname    : File name. If defined, will save the figure and close it.
      cmap     : Colormap to use (2D data only). See matplotlib doc for
                 options.
      labels   : Labels to add to the figure. {'label':position (deg)}
                 dictionnary.
      log      : Logaritmic color scale (def. False)
      vmin     : Minimal value of the colorscale (2D data only).
      vmax     : Maximal value of the colorscale (2D data only).
      showpts  : If set to True, will show datapoints.
    """
    if "showpts" in kwargs and kwargs["showpts"]:
        xv, yv = _np.meshgrid(_np.deg2rad(phi), 90 - _np.array(r))
    if "interp" in kwargs and kwargs["interp"] != "None":
        ndata = _np.concatenate([data, data[:, 0, ...][:, None]], axis=1)
        nphi = _np.linspace(phi[0], phi[0] + 360, n * len(phi), endpoint=False)
        nr = _np.linspace(r[0], r[-1], n * (len(r) - 1) + 1)
        if data.ndim == 3:
            data = _np.zeros((len(nr), len(nphi), data.shape[2]))
            for i in range(data.shape[2]):
                interp = _interp2d(
                    _np.concatenate([phi, [phi[0] + 360]]),
                    r,
                    ndata[:, :, i],
                    kind=kwargs["interp"],
                )
                data[:, :, i] = interp(nphi, nr)
        else:
            interp = _interp2d(
                _np.concatenate([phi, [phi[0] + 360]]), r, ndata, kind=kwargs["interp"]
            )
            data = interp(nphi, nr)
        phi = nphi
        r = nr
        n = 1
    else:
        data = _np.repeat(data, n, axis=1)

    r = _np.asarray(r).tolist()
    r[0:0] = [2 * r[0] - r[1]]
    r[-1:-1] = [r[-1]]
    r = 90 - _np.mean([r[:-1], r[1:]], 0)
    r[r < 0] = 0
    r[r > 90] = 90

    theta = _np.linspace(0, 2 * _np.pi, n * len(phi) + 1) - _np.mean(
        _np.radians(phi[:2])
    )

    Theta, R = _np.meshgrid(theta, r)
    if "autogain" in kwargs and kwargs["autogain"]:
        gain = _np.max(data)
        data /= gain
    if data.ndim == 3:
        color = data.reshape((-1, data.shape[2]))

    _plt.figure()
    ax = _plt.subplot(111, polar=True)
    ax.set_theta_zero_location("N")
    ax.xaxis.set_ticklabels(["N", "NE", "E", "SE", "S", "SW", "W", "NW"])
    if data.ndim == 3:
        m = _plt.pcolormesh(
            Theta, R, data[:, :, 0], color=color, linewidth=0, vmin=0, vmax=1
        )
        m.set_array(None)
    else:
        args = dict()
        if "cmap" in kwargs:
            args["cmap"] = kwargs["cmap"]
        if "vmin" in kwargs:
            args["vmin"] = kwargs["vmin"]
        if "vmax" in kwargs:
            args["vmax"] = kwargs["vmax"]
        if "log" in kwargs and kwargs["log"]:
            args["norm"] = _mpl.colors.LogNorm(vmin=kwargs["vmin"], vmax=kwargs["vmax"])
        m = _plt.pcolormesh(Theta, R, data, linewidth=0, **args)

    if "showpts" in kwargs and kwargs["showpts"]:
        _plt.plot(xv, yv, "r.")
    title_str = ""
    if "title" in kwargs:
        title_str += kwargs["title"]
    if "autogain" in kwargs and kwargs["autogain"]:
        if title_str != "":
            title_str += "\n"
        title_str += "gain = %.4g" % gain
    if title_str != "":
        _plt.title(title_str)

    if "clabel" in kwargs:
        _plt.colorbar(label=kwargs["clabel"], pad=0.11)

    if "labels" in kwargs:
        for name, pos in list(kwargs["labels"].items()):
            _plt.annotate(
                name,
                xy=[_np.radians(pos), _np.max(r)],
                xytext=[_np.radians(pos), _np.max(r) + 7],
                verticalalignment="top" if (90 < pos < 270) else "bottom",
                horizontalalignment="right" if pos < 180 else "left",
                annotation_clip=False,
                arrowprops=dict(facecolor="black", width=0.1, headlength=0.1),
            )

    _plt.tight_layout()

    if "fname" in kwargs:
        _plt.savefig(kwargs["fname"])
        _plt.close()
