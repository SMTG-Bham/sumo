# coding: utf-8
# Copyright (c) Scanlon Materials Theory Group
# Distributed under the terms of the MIT License.

"""
Subpackage providing helper functions for generating publication ready plots.
"""

import numpy as np

import matplotlib.pyplot
from matplotlib.collections import LineCollection
from matplotlib import rc, rcParams

colour_cache = {}


def pretty_plot(width=None, height=None, plt=None, dpi=None, fonts=None):
    """Get a :obj:`matplotlib.pyplot` object with publication ready defaults.

    Args:
        width (:obj:`float`, optional): The width of the plot.
        height (:obj:`float`, optional): The height of the plot.
        plt (:obj:`matplotlib.pyplot`, optional): A :obj:`matplotlib.pyplot`
            object to use for plotting.
        dpi (:obj:`int`, optional): The dots-per-inch (pixel density) for
            the plot.
        fonts (:obj:`list`, optional): Fonts to use in the plot. Can be a
            a single font, specified as a :obj:`str`, or several fonts,
            specified as a :obj:`list` of :obj:`str`.

    Returns:
        :obj:`matplotlib.pyplot`: A :obj:`matplotlib.pyplot` object with
        publication ready defaults set.
    """

    if plt is None:
        plt = matplotlib.pyplot
        if width is None:
            width = matplotlib.rcParams['figure.figsize'][0]
        if height is None:
            height = matplotlib.rcParams['figure.figsize'][1]

        if dpi is not None:
            matplotlib.rcParams['figure.dpi'] = dpi

        fig = plt.figure(figsize=(width, height))
        fig.add_subplot(1, 1, 1)

    if fonts is not None:
        if type(fonts) is str:
            fonts = [fonts]

        rc('font', **{'family': 'sans-serif', 'sans-serif': fonts})
        rc('text', usetex=False)

    return plt


def pretty_subplot(nrows, ncols, width=None, height=None, sharex=True,
                   sharey=True, dpi=None, fonts=None, plt=None,
                   gridspec_kw=None):
    """Get a :obj:`matplotlib.pyplot` subplot object with pretty defaults.

    Args:
        nrows (int): The number of rows in the subplot.
        ncols (int): The number of columns in the subplot.
        width (:obj:`float`, optional): The width of the plot.
        height (:obj:`float`, optional): The height of the plot.
        sharex (:obj:`bool`, optional): All subplots share the same x-axis.
            Defaults to ``True``.
        sharey (:obj:`bool`, optional): All subplots share the same y-axis.
            Defaults to ``True``.
        dpi (:obj:`int`, optional): The dots-per-inch (pixel density) for
            the plot.
        plt (:obj:`matplotlib.pyplot`, optional): A :obj:`matplotlib.pyplot`
            object to use for plotting.
        fonts (:obj:`list`, optional): Fonts to use in the plot. Can be a
            a single font, specified as a :obj:`str`, or several fonts,
            specified as a :obj:`list` of :obj:`str`.
        gridspec_kw (:obj:`dict`, optional): Gridspec parameters. Please see:
            :obj:`matplotlib.pyplot.subplot` for more information. Defaults
            to ``None``.

    Returns:
        :obj:`matplotlib.pyplot`: A :obj:`matplotlib.pyplot` subplot object
        with publication ready defaults set.
    """

    if width is None:
        width = rcParams['figure.figsize'][0]
    if height is None:
        height = rcParams['figure.figsize'][1]

    # TODO: Make this work if plt is already set...
    if plt is None:
        plt = matplotlib.pyplot
        f, axes = plt.subplots(nrows, ncols, sharex=sharex, sharey=sharey,
                               dpi=dpi, figsize=(width, height), facecolor='w',
                               gridspec_kw=gridspec_kw)

    if fonts is not None:
        if type(fonts) is str:
            fonts = [fonts]

        rc('font', **{'family': 'sans-serif', 'sans-serif': fonts})
        rc('text', usetex=False)

    rc('legend', handlelength=1.5)
    return plt


def power_tick(val, pos):
    """Custom power ticker function. """
    if val == 0:
        return r'$\mathregular{0}$'
    exponent = int(np.log10(val))
    coeff = val / 10**exponent
    return r'$\mathregular{{{:0.1f} x 10^{:2d}}}$'.format(coeff, exponent)


def rgbline(x, y, red, green, blue, alpha=1, linestyles="solid",
            linewidth=2.5):
    """Get a RGB coloured line for plotting.

    Args:
        x (list): x-axis data.
        y (list): y-axis data (can be multidimensional array).
        red (list): Red data (must have same shape as ``y``).
        green (list): Green data (must have same shape as ``y``).
        blue (list): blue data (must have same shape as ``y``).
        alpha (:obj:`list` or :obj:`int`, optional): Alpha (transparency)
            data (must have same shape as ``y`` or be an :obj:`int`).
        linestyles (:obj:`str`, optional): Linestyle for plot. Options are
            ``"solid"`` or ``"dotted"``.
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
