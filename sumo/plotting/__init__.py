# Copyright (c) Scanlon Materials Theory Group
# Distributed under the terms of the MIT License.

"""
Subpackage providing helper functions for generating publication ready plots.
"""
import os
from functools import wraps

import matplotlib.pyplot
import numpy as np
from matplotlib import rcParams
from matplotlib.collections import LineCollection

try:
    from importlib.resources import files as ilr_files
except ImportError:  # Python < 3.9
    from importlib_resources import files as ilr_files

colour_cache = {}

sumo_base_style = ilr_files("sumo.plotting") / "sumo_base.mplstyle"
sumo_dos_style = ilr_files("sumo.plotting") / "sumo_dos.mplstyle"
sumo_bs_style = ilr_files("sumo.plotting") / "sumo_bs.mplstyle"
sumo_phonon_style = ilr_files("sumo.plotting") / "sumo_phonon.mplstyle"
sumo_optics_style = ilr_files("sumo.plotting") / "sumo_optics.mplstyle"


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
        @wraps(get_plot)
        def wrapper(*args, fonts=None, style=None, no_base_style=False, **kwargs):
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
                list_style += [{"font.family": "sans-serif", "font.sans-serif": fonts}]

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
            width = matplotlib.rcParams["figure.figsize"][0]
        if height is None:
            height = matplotlib.rcParams["figure.figsize"][1]

        if dpi is not None:
            matplotlib.rcParams["figure.dpi"] = dpi

        fig = plt.figure(figsize=(width, height))
        fig.add_subplot(1, 1, 1)

    return plt


def pretty_subplot(
    nrows,
    ncols,
    width=None,
    height=None,
    sharex=True,
    sharey=True,
    dpi=None,
    plt=None,
    gridspec_kw=None,
):
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
        width = rcParams["figure.figsize"][0]
    if height is None:
        height = rcParams["figure.figsize"][1]

    # TODO: Make this work if plt is already set...
    if plt is None:
        plt = matplotlib.pyplot
        plt.subplots(
            nrows,
            ncols,
            sharex=sharex,
            sharey=sharey,
            dpi=dpi,
            figsize=(width, height),
            facecolor="w",
            gridspec_kw=gridspec_kw,
        )

    return plt


def curry_power_tick(times_sign=r"\times"):
    def f(val, pos):
        return power_tick(val, pos, times_sign=times_sign)

    return f


def power_tick(val, pos, times_sign=r"\times"):
    """Custom power ticker function."""
    if val == 0:
        return r"$\mathregular{0}$"
    elif val < 0:
        exponent = int(np.log10(-val))
    else:
        exponent = int(np.log10(val))
    coeff = val / 10**exponent
    prec = 0 if coeff % 1 == 0 else 1

    return rf"${coeff:.{prec}f}\mathrm{{{times_sign}}}10^{{{exponent:2d}}}$"


def colorline(
    x,
    y,
    weights,
    color1="#FF0000",
    color2="#00FF00",
    color3="#0000FF",
    colorspace="lab",
    linestyles="solid",
    linewidth=2.5,
):
    """Get a RGB coloured line for plotting.

    Args:
        x (list): x-axis data.
        y (list): y-axis data (can be multidimensional array).
        weights (list): The weights of the color1, color2, and color3 channels.
            Given as an array with the shape (n, 3), where n is the same length
            as the x and y data.
        color1 (str): A color specified in any way supported by matplotlib.
        color2 (str): A color specified in any way supported by matplotlib.
        color3 (str): A color specified in any way supported by matplotlib.
        colorspace (str): The colorspace in which to perform the interpolation.
            The allowed values are rgb, hsv, lab, luvlc, lablch, and xyz.
        linestyles (:obj:`str`, optional): Linestyle for plot. Options are
            ``"solid"`` or ``"dotted"``.
    """
    y = np.array(y)
    if len(y.shape) == 1:
        y = np.array([y])
        weights = np.array([weights])

    seg = []
    colours = []
    for yy, ww in zip(y, weights):
        pts = np.array([x, yy]).T.reshape(-1, 1, 2)
        if len(pts) > 1:  # need at least one point to interpolate colours
            seg.extend(np.concatenate([pts[:-1], pts[1:]], axis=1))

            nseg = len(x) - 1
            w = [0.5 * (ww[i] + ww[i + 1]) for i in range(nseg)]
            c = get_interpolated_colors(
                color1, color2, color3, w, colorspace=colorspace
            )
            colours.extend(c.tolist())

    lc = LineCollection(
        seg,
        colors=colours,
        rasterized=True,
        linewidth=linewidth,
        linestyles=linestyles,
    )
    return lc


def get_interpolated_colors(color1, color2, color3, weights, colorspace="lab"):
    """
    Interpolate colors at a number of points within a colorspace.

    Args:
        color1 (str): A color specified in any way supported by matplotlib.
        color2 (str): A color specified in any way supported by matplotlib.
        color3 (str): A color specified in any way supported by matplotlib.
        weights (list): A list of weights with the shape (n, 3).
            Where the 3 values of the last axis give the amount of
            color1, color2, and color3.
        colorspace (str): The colorspace in which to perform the interpolation.
            The allowed values are rgb, hsv, lab, luvlc, lablch, and xyz.

    Returns:
        A list of colors, specified in the rgb format as a (n, 3) array.
    """
    from colormath.color_conversions import convert_color
    from colormath.color_objects import (
        HSVColor,
        LabColor,
        LCHabColor,
        LCHuvColor,
        XYZColor,
        sRGBColor,
    )
    from matplotlib.colors import to_rgb

    colorspace_mapping = {
        "rgb": sRGBColor,
        "hsv": HSVColor,
        "lab": LabColor,
        "luvlch": LCHuvColor,
        "lablch": LCHabColor,
        "xyz": XYZColor,
    }
    if colorspace not in list(colorspace_mapping.keys()):
        raise ValueError(f"colorspace must be one of {colorspace_mapping.keys()}")

    colorspace = colorspace_mapping[colorspace]

    # first convert matplotlib color specification to colormath sRGB
    color1_rgb = sRGBColor(*to_rgb(color1))
    color2_rgb = sRGBColor(*to_rgb(color2))
    color3_rgb = sRGBColor(*to_rgb(color3))

    # now convert to the colorspace basis for interpolation
    basis1 = np.array(
        convert_color(color1_rgb, colorspace, target_illuminant="d50").get_value_tuple()
    )
    basis2 = np.array(
        convert_color(color2_rgb, colorspace, target_illuminant="d50").get_value_tuple()
    )
    basis3 = np.array(
        convert_color(color3_rgb, colorspace, target_illuminant="d50").get_value_tuple()
    )

    # ensure weights is a numpy array
    weights = np.asarray(weights)

    # perform the interpolation in the colorspace basis
    colors = (
        basis1 * weights[:, 0][:, None]
        + basis2 * weights[:, 1][:, None]
        + basis3 * weights[:, 2][:, None]
    )

    # convert colors to RGB
    rgb_colors = [
        convert_color(colorspace(*c), sRGBColor).get_value_tuple() for c in colors
    ]

    # ensure all rgb values are less than 1 (sometimes issues in interpolation
    # gives values slightly over 1)
    return np.minimum(rgb_colors, 1)


def draw_themed_line(y, ax, orientation="horizontal", **kwargs):
    """Draw a horizontal line using the theme settings

    Args:
        y (float): Position of line in data coordinates
        ax (Axes): Matplotlib Axes on which line is drawn
        orientation (str, optional): Orientation of line. Options are
            ``"horizontal"`` or ``"vertical"``.
        **kwargs: Additional keyword arguments passed to ``ax.axhline`` or
            ``ax.axvline``, which can be used to override the theme settings.
    """

    # Note to future developers: feel free to add plenty more optional
    # arguments to this to mess with linestyle, zorder etc.
    # Just .update() the options dict

    themed_line_options = dict(
        color=rcParams["grid.color"],
        linestyle="--",
        dashes=(5, 2),
        zorder=0,
        linewidth=rcParams["ytick.major.width"],
    )
    themed_line_options.update(kwargs)

    if orientation == "horizontal":
        ax.axhline(y, **themed_line_options)
    elif orientation == "vertical":
        ax.axvline(y, **themed_line_options)
    else:
        raise ValueError(f'Line orientation "{orientation}" not supported')
