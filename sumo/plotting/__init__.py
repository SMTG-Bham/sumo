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
from pkg_resources import resource_filename

colour_cache = {}

sumo_base_style = resource_filename('sumo.plotting', 'sumo_base.mplstyle')
sumo_dos_style = resource_filename('sumo.plotting', 'sumo_dos.mplstyle')
sumo_bs_style = resource_filename('sumo.plotting', 'sumo_bs.mplstyle')
sumo_phonon_style = resource_filename('sumo.plotting', 'sumo_phonon.mplstyle')
sumo_optics_style = resource_filename('sumo.plotting', 'sumo_optics.mplstyle')


def styled_plot(*style_sheets):
    """Return a decorator that will apply matplotlib style sheets to a plot.

    ``style_sheets`` is a base set of styles, which will be ignored if
    ``no_base_style`` is set in the decorated function arguments.

    The style will further be overwritten by any styles in the ``style``
    optional argument of the decorated function.

    Args:
        style_sheets (:obj:`list`, :obj:`str`, or :obj:`dict`): Any matplotlib
            supported definition of a style sheet. Can be a list of style of
            style sheets.
    """

    def decorator(get_plot):

        def wrapper(*args, fonts=None, style=None, no_base_style=False,
                    **kwargs):

            if no_base_style:
                list_style = []
            else:
                list_style = list(style_sheets)

            if style is not None:
                if isinstance(style, list):
                    list_style += style
                else:
                    list_style += [style]

            if fonts is not None:
                list_style += [{'font.family': 'sans-serif',
                               'font.sans-serif': fonts}]

            matplotlib.pyplot.style.use(list_style)
            return get_plot(*args, **kwargs)

        return wrapper
    return decorator


def pretty_plot(width=None, height=None, plt=None, dpi=None):
    """Get a :obj:`matplotlib.pyplot` object with publication ready defaults.

    Args:
        width (:obj:`float`, optional): The width of the plot.
        height (:obj:`float`, optional): The height of the plot.
        plt (:obj:`matplotlib.pyplot`, optional): A :obj:`matplotlib.pyplot`
            object to use for plotting.
        dpi (:obj:`int`, optional): The dots-per-inch (pixel density) for
            the plot.

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

    return plt


def pretty_subplot(nrows, ncols, width=None, height=None, sharex=True,
                   sharey=True, dpi=None, plt=None, gridspec_kw=None):
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
        plt.subplots(nrows, ncols, sharex=sharex, sharey=sharey, dpi=dpi,
                     figsize=(width, height), facecolor='w',
                     gridspec_kw=gridspec_kw)

    return plt


def curry_power_tick(times_sign=r'\times'):
    def f(val, pos):
        return power_tick(val, pos, times_sign=times_sign)
    return f

def power_tick(val, pos, times_sign=r'\times'):
    """Custom power ticker function. """
    if val == 0:
        return r'$\mathregular{0}$'
    elif val < 0:
        exponent = int(np.log10(-val))
    else:
        exponent = int(np.log10(val))
    coeff = val / 10**exponent

    return r'$\mathregular{{{:.1f} {} 10^{:2d}}}$'.format(coeff,
                                                          times_sign,
                                                          exponent)

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
    elif isinstance(alpha, int):
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
