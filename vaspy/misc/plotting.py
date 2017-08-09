#! /usr/bin/env python

from itertools import cycle

default_colours = [[240, 163, 255], [0, 117, 220], [153, 63, 0], [76, 0, 92],
                   [66, 102, 0], [255, 0, 16], [157, 204, 0], [194, 0, 136],
                   [0, 51, 128], [255, 164, 5], [255, 255, 0], [255, 80, 5],
                   [94, 241, 242], [116, 10, 255], [153, 0, 0], [0, 153, 143],
                   [0, 92, 49], [43, 206, 72], [255, 204, 153], [148, 255, 181],
                   [143, 124, 0], [255, 168, 187], [128, 128, 128]]

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
                   labelsize=_ticklabelsize, pad=10, direction='in')
    ax.tick_params(which='minor', size=_ticksize/2, width=_linewidth,
                   direction='in')

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
                       labelsize=_ticklabelsize, pad=10, direction='in')
        ax.tick_params(which='minor', size=_ticksize/2, width=_linewidth,
                       direction='in')

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
    rc('legend', handlelength=1)
    return plt


def colour_cycle():
    from numpy import array
    rgb_colours = array(default_colours)/255.
    return cycle(rgb_colours)


def colour_cycler():
    from numpy import array
    from cycler import cycler
    rgb_colours = array(default_colours)/255.
    return cycler('color', rgb_colours)
