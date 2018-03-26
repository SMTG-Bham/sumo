# coding: utf-8
# Copyright (c) Scanlon Materials Theory Group
# Distributed under the terms of the MIT License.

import numpy as np

from cycler import cycler
from itertools import cycle
from matplotlib.collections import LineCollection

default_colours = [[240, 163, 255], [0, 117, 220], [153, 63, 0], [76, 0, 92],
                   [66, 102, 0], [255, 0, 16], [157, 204, 0], [194, 0, 136],
                   [0, 51, 128], [255, 164, 5], [255, 255, 0], [255, 80, 5],
                   [94, 241, 242], [116, 10, 255], [153, 0, 0], [0, 153, 143],
                   [0, 92, 49], [43, 206, 72], [255, 204, 153],
                   [148, 255, 181], [143, 124, 0], [255, 168, 187],
                   [128, 128, 128]]

default_fonts = ['Whitney Book Extended', 'Arial', 'Whitney Book', 'Helvetica',
                 'Liberation Sans', 'Andale Sans']

_ticklabelsize = 22
_labelsize = 22
_ticksize = 15
_linewidth = 1.3


def pretty_plot(width=5, height=5, plt=None, dpi=None, fonts=None):
    from matplotlib import rc

    if plt is None:
        import matplotlib.pyplot as plt
        plt.figure(figsize=(width, height), facecolor="w", dpi=dpi)
        ax = plt.gca()
        ax.set_prop_cycle(colour_cycler())

    ax = plt.gca()

    ax.tick_params(width=_linewidth, size=_ticksize)
    ax.tick_params(which='major', size=_ticksize, width=_linewidth,
                   labelsize=_ticklabelsize, pad=7, direction='in',
                   right='off', top='off')
    ax.tick_params(which='minor', size=_ticksize/2, width=_linewidth,
                   direction='in', right='off', top='off')

    ax.set_title(ax.get_title(), size=20)
    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(_linewidth)

    ax.set_xlabel(ax.get_xlabel(), size=_labelsize)
    ax.set_ylabel(ax.get_ylabel(), size=_labelsize)

    fonts = default_fonts if fonts is None else fonts + default_fonts

    rc('font', **{'family': 'sans-serif', 'sans-serif': fonts})
    rc('text', usetex=False)
    rc('pdf', fonttype=42)
    rc('mathtext', fontset='stixsans')
    rc('legend', handlelength=2)
    return plt


def pretty_subplot(nrows, ncols, width=5, height=5, sharex=True,
                   sharey=True, dpi=None, fonts=None, plt=None,
                   gridspec_kw=None):
    from matplotlib import rc

    # TODO: Make this work if plt is already set...
    if plt is None:
        import matplotlib.pyplot as plt
        f, axes = plt.subplots(nrows, ncols, sharex=sharex, sharey=sharey,
                               dpi=dpi, figsize=(width, height), facecolor='w',
                               gridspec_kw=gridspec_kw)

    for ax in axes:
        ax.set_prop_cycle(colour_cycler())
        ax.tick_params(width=_linewidth, size=_ticksize)
        ax.tick_params(which='major', size=_ticksize, width=_linewidth,
                       labelsize=_ticklabelsize, pad=7, direction='in',
                       right='off', top='off')
        ax.tick_params(which='minor', size=_ticksize/2, width=_linewidth,
                       direction='in', right='off', top='off')

        ax.set_title(ax.get_title(), size=20)
        for axis in ['top', 'bottom', 'left', 'right']:
            ax.spines[axis].set_linewidth(_linewidth)

        ax.set_xlabel(ax.get_xlabel(), size=_labelsize)
        ax.set_ylabel(ax.get_ylabel(), size=_labelsize)

    fonts = default_fonts if fonts is None else fonts + default_fonts

    rc('font', **{'family': 'sans-serif', 'sans-serif': fonts})
    rc('text', usetex=False)
    rc('pdf', fonttype=42)
    rc('mathtext', fontset='stixsans')
    rc('legend', handlelength=1.5)
    return plt


def colour_cycle():
    rgb_colours = np.array(default_colours)/255.
    return cycle(rgb_colours)


def colour_cycler():
    rgb_colours = np.array(default_colours)/255.
    return cycler('color', rgb_colours)


def power_tick(val, pos):
    if val == 0:
        return '$\mathregular{0}$'
    exponent = int(np.log10(val))
    coeff = val / 10**exponent
    return '$\mathregular{{{:0.1f} x 10^{:2d}}}$'.format(coeff, exponent)


def rgbline(x, y, red, green, blue, alpha=1, linestyles="solid", linewidth=2.5):
    """An RGB colored line for plotting.

    Args:
        ax: matplotlib axis
        x: x-axis data
        y: y-axis data (can be multidimensional array)
        red: red data (must have same shape as y)
        green: green data (must have same shape as y)
        blue: blue data (must have same shape as y)
        alpha: alpha values data (must have same shape as y or be int)
        linestyles: linestyle for plot (e.g., "solid" or "dotted")
    """
    y = np.array(y)
    if len(y.shape) == 1:
        y = np.array([y])
        red = np.array([red])
        green = np.array([green])
        blue = np.array([blue])
        alpha = np.array([alpha])
    elif type(alpha) == int:
        alpha = [alpha] * len(y)

    seg = []
    colours = []
    for yy, rr, gg, bb, aa in zip(y, red, green, blue, alpha):
        pts = np.array([x, yy]).T.reshape(-1, 1, 2)
        seg.extend(np.concatenate([pts[:-1], pts[1:]], axis=1))

        nseg = len(x) - 1
        r = [0.5 * (rr[i] + rr[i + 1]) for i in range(nseg)]
        g = [0.5 * (gg[i] + gg[i + 1]) for i in range(nseg)]
        b = [0.5 * (bb[i] + bb[i + 1]) for i in range(nseg)]
        a = np.ones(nseg, np.float) * aa
        colours.extend(list(zip(r, g, b, a)))

    lc = LineCollection(seg, colors=colours, rasterized=True,
                        linewidth=linewidth, linestyles=linestyles)
    return lc
