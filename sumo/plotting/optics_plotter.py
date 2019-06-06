# coding: utf-8
# Copyright (c) Scanlon Materials Theory Group
# Distributed under the terms of the MIT License.

"""
This module provides a class for plotting optical absorption spectra.
"""

import numpy as np

from matplotlib import rcParams
from matplotlib.ticker import MaxNLocator, FuncFormatter, AutoMinorLocator
from matplotlib.font_manager import findfont, FontProperties

from sumo.plotting import (pretty_subplot, curry_power_tick, styled_plot,
                           sumo_base_style, sumo_optics_style)


class SOpticsPlotter(object):
    """Class for plotting optical spectra.

    The easiest way to initialise this class is using the output from the
    :obj:`sumo.electronic_structure.optics.calculate_alpha()` method.

    Args:
        spec_data (:obj:`dict` or :obj:`tuple` or :obj:`list`):
            The optical absorption spectra. Should be formatted as a dict of
            :obj:`tuple` of :obj:`list`::

                {property: ([energies], [alpha]), ...}

            Alternatively, the anisotropic (directional dependent) absorption
            can be plotted if the data formatted as::

                {property: ([energies], [alpha_xx, alpha_yy, alpha_zz]), ...}

            If a :obj:`list` of :obj:`tuple` is provided, then multiple
            absorption spectra can be plotted simultaneously.

            The *property* keys set the type of spectrum being plotted (and
            determine the y-axis label). Recognised values are 'absorption',
            'loss', 'eps-real', 'eps-im'. Other values will be plotted and used
            as the y-axis label.

            It is recommended to use an OrderedDict to get a predictable
            sequence of plotting axes (from top to bottom).

            If no dict at all is used (i.e. spec_data is a data tuple or a list
            of tuples) property will default to 'abs'

        band_gap (:obj:`float` or :obj:`list`, optional): The band gap as a
            :obj:`float`, plotted as a dashed line. If plotting multiple
            spectra then a :obj:`list` of band gaps can be provided. If an item
            in the list is :obj:`None` no corresponding line is plotted.
        label (:obj:`str` or :obj:`list`): A label to identify the spectra.
            If plotting multiple spectra then a :obj:`list` of labels can
            be provided.
    """

    def __init__(self, spec_data, band_gap=None, label=None):
        # Standardise to dict/list/tuple format
        if isinstance(spec_data, tuple):
            spec_data = {'absorption': [spec_data]}
        elif isinstance(spec_data, list):
            spec_data = {'absorption': spec_data}
        for k, v in spec_data.items():
            spec_data[k] = (v if isinstance(v, list) else [v])

        n_datasets = len(next(iter(spec_data.items()))[1])

        if isinstance(band_gap, float):
            band_gap = [band_gap]
        elif not band_gap:
            band_gap = [None] * n_datasets

        if isinstance(label, str):
            label = [label]
        elif not label and n_datasets > 1:
            label = [str(i) for i in range(1, n_datasets + 1)]
        elif not label:
            label = [''] * n_datasets

        if len(set([len(label), n_datasets, len(band_gap)])) > 1:
            raise ValueError('spec_data, band_gap, and label not same size')

        if 'absorption' in spec_data:
            # xmax: find lowest energy where absorption > 2e4; add 1 eV
            xmax = 0
            for ener, alpha in spec_data['absorption']:
                    mask = alpha > 2e4
                    if len(alpha.shape) == 1:
                        x = min(ener[mask])
                    else:
                        x = max([min(ener[mask[:, i]]) for i in range(3)])
                    xmax = x if x > xmax else xmax
        else:
            # If no absorption data, use the maximum of the shortest x range
            xmax = float('Inf')
            for _, spec in spec_data.items():
                for dataset in spec:
                    xmax = min((xmax, np.max(dataset[0])))

        self._spec_data = spec_data
        self._band_gap = band_gap
        self._label = label
        self._xmax = xmax + 1.

    @styled_plot(sumo_base_style, sumo_optics_style)
    def get_plot(self, width=None, height=None, xmin=0., xmax=None, ymin=0,
                 ymax=None, colours=None, dpi=400, plt=None, fonts=None,
                 style=None, no_base_style=False):
        """Get a :obj:`matplotlib.pyplot` object of the optical spectra.

        Args:
            width (:obj:`float`, optional): The width of the plot.
            height (:obj:`float`, optional): The height of the plot.
            xmin (:obj:`float`, optional): The minimum energy on the x-axis.
            xmax (:obj:`float`, optional): The maximum energy on the x-axis.
            ymin (:obj:`float`, optional): The minimum absorption intensity on
                the y-axis.
            ymax (:obj:`float`, optional): The maximum absorption intensity on
                the y-axis.
            colours (:obj:`list`, optional): A :obj:`list` of colours to use in
                the plot. The colours can be specified as a hex code, set of
                rgb values, or any other format supported by matplotlib.
            dpi (:obj:`int`, optional): The dots-per-inch (pixel density) for
                the image.
            plt (:obj:`matplotlib.pyplot`, optional): A
                :obj:`matplotlib.pyplot` object to use for plotting.
            fonts (:obj:`list`, optional): Fonts to use in the plot. Can be a
                a single font, specified as a :obj:`str`, or several fonts,
                specified as a :obj:`list` of :obj:`str`.
            style (:obj:`list`, :obj:`str`, or :obj:`dict`): Any matplotlib
                style specifications, to be composed on top of Sumo base
                style.
            no_base_style (:obj:`bool`, optional): Prevent use of sumo base
                style. This can make alternative styles behave more
                predictably.

        Returns:
            :obj:`matplotlib.pyplot`: The plot of optical spectra.
        """
        n_plots = len(self._spec_data)
        plt = pretty_subplot(n_plots, 1, sharex=True, sharey=False,
                             width=width, height=height,
                             dpi=dpi, plt=plt)
        fig = plt.gcf()

        optics_colours = rcParams['axes.prop_cycle'].by_key()['color']
        if colours is not None:
            optics_colours = colours + optics_colours

        standard_ylabels = {
            'absorption': r'Absorption (cm$^\mathregular{-1}$)',
            'loss': r'Energy-loss',
            'eps_real': r'Re($\epsilon$)',
            'eps_imag': r'Im($\epsilon$)',
            'n_real': r'Re(n)',
            'n_imag': r'Im(n)'}

        if ymax is None:
            ymax_series = [None] * n_plots
        elif isinstance(ymax, float) or isinstance(ymax, int):
            ymax_series = [ymax] * n_plots
        elif not isinstance(ymax, list):
            raise ValueError()
        else:
            ymax_series = ymax

        if ymin is None:
            ymin_series = [None] * n_plots
        elif isinstance(ymin, float) or isinstance(ymin, int):
            ymin_series = [ymin] * n_plots
        elif not isinstance(ymin, list):
            raise ValueError()
        else:
            ymin_series = ymin

        for i, (spectrum_key, data), ymin, ymax in zip(range(n_plots),
                                                       self._spec_data.items(),
                                                       ymin_series,
                                                       ymax_series):
            ax = fig.axes[i]
            _plot_spectrum(data, self._label, self._band_gap,
                           ax, optics_colours)

            xmax = xmax if xmax else self._xmax
            ax.set_xlim(xmin, xmax)

            if ymin is None and spectrum_key in ('absorption', 'loss',
                                                 'eps_imag', 'n_imag'):
                ymin = 0
            elif ymin is None:
                ymin = ax.get_ylim()[0]

            if ymax is None and spectrum_key in ('absorption',):
                ymax = 1e5
            elif ymax is None:
                ymax = ax.get_ylim()[1]

            ax.set_ylim(ymin, ymax)

            if spectrum_key == 'absorption':
                font = findfont(FontProperties(family=['sans-serif']))
                if 'Whitney' in font:
                    times_sign = 'x'
                else:
                    times_sign = r'\times'
                ax.yaxis.set_major_formatter(
                    FuncFormatter(curry_power_tick(times_sign=times_sign)))

            ax.yaxis.set_major_locator(MaxNLocator(5))
            ax.xaxis.set_major_locator(MaxNLocator(3))
            ax.yaxis.set_minor_locator(AutoMinorLocator(2))
            ax.xaxis.set_minor_locator(AutoMinorLocator(2))

            ax.set_ylabel(standard_ylabels.get(spectrum_key, spectrum_key))

            if i == 0:
                if (not np.all(np.array(self._label) == '') or
                        len(np.array(next(iter(
                            self._spec_data.items()))[1][0][1]).shape) > 1):
                    ax.legend(loc='best', frameon=False, ncol=1)

        ax.set_xlabel('Energy (eV)')

        # If only one plot, fix aspect ratio to match canvas
        if len(self._spec_data) == 1:
            x0, x1 = ax.get_xlim()
            y0, y1 = ax.get_ylim()
            if width is None:
                width = rcParams['figure.figsize'][0]
            if height is None:
                height = rcParams['figure.figsize'][1]
            ax.set_aspect((height/width) * ((x1-x0)/(y1-y0)))

        # Otherwise, rely only on tight_layout and hope for the best
        plt.tight_layout()
        return plt


def _plot_spectrum(data, label, band_gap, ax, optics_colours):
    for (ener, alpha), abs_label, bg, c in zip(data,
                                               label,
                                               band_gap,
                                               optics_colours):
        if len(alpha.shape) == 1:
            # if averaged optics only plot one line
            ax.plot(ener, alpha, label=abs_label, c=c)

        else:
            data = zip(range(3), ['xx', 'yy', 'zz'], ['-', '--', '-.'])

            for direction_id, direction_label, ls in data:
                if not abs_label:
                    label = direction_label
                else:
                    label = r'{}$_\mathregular{{{}}}$'
                    label = label.format(abs_label, direction_label)

                ax.plot(ener, alpha[:, direction_id], ls=ls,
                        label=label, c=c)
        if bg:
            # plot band gap line
            ax.axvline(bg, ls=':', c=c)
