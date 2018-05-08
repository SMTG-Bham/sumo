# coding: utf-8
# Copyright (c) Scanlon Materials Theory Group
# Distributed under the terms of the MIT License.

"""
This module provides a class for plotting density of states data.
"""

import itertools

from sumo.electronic_structure.dos import sort_orbitals
from sumo.plotting import pretty_plot, pretty_subplot, colour_cycle
from sumo.plotting import colour_cache

from pymatgen.electronic_structure.core import Spin

try:
    import configparser
except ImportError:
    import ConfigParser as configparser

line_width = 1.5
empty_space = 1.05
label_size = 22
band_linewidth = 2
col_cycle = colour_cycle()


class SDOSPlotter(object):
    """Class for plotting density of states data.

    This class should be initialised with the total and partial density
    of states. The easiest way to generate the partial density of states is
    using the following method::

        pdos = sumo.electronic_structure.dos.get_pdos()

    Args:
        dos (:obj:`~pymatgen.electronic_structure.dos.Dos`): The total density
            of states.
        pdos (:obj:`dict`, optional): The partial density of states. Formatted
            as a :obj:`dict` of :obj:`dict` mapping the elements and their
            orbitals to :obj:`~pymatgen.electronic_structure.dos.Dos` objects.
            For example::

                {
                    'Bi': {'s': Dos, 'p': Dos ... },
                    'S': {'s': Dos}
                }

            Usually generated using the
            :obj:`sumo.electronric_structure.dos.get_pdos()` function.
    """

    def __init__(self, dos, pdos=None):
        self._dos = dos
        self._pdos = pdos

    def dos_plot_data(self, yscale=1, xmin=-6., xmax=6., colours=None,
                      plot_total=True, legend_cutoff=3, subplot=False):
        """Get the plotting data.

        Args:
            yscale (:obj:`float`, optional): Scaling factor for the y-axis.
            xmin (:obj:`float`, optional): The minimum energy to mask the
                energy and density of states data (reduces plotting load).
            xmax (:obj:`float`, optional): The maximum energy to mask the
                energy and density of states data (reduces plotting load).
            colours (:obj:`dict`, optional): Use custom colours for specific
                element and orbital combinations. Specified as a :obj:`dict` of
                :obj:`dict` of the colours. For example::

                    {
                        'Sn': {'s': 'r', 'p': 'b'},
                        'O': {'s': '#000000'}
                    }

                The colour can be a hex code, series of rgb value, or any other
                format supported by matplotlib.
            plot_total (:obj:`bool`, optional): Plot the total density of
                states. Defaults to ``True``.
            legend_cutoff (:obj:`float`, optional): The cut-off (in % of the
                maximum density of states within the plotting range) for an
                elemental orbital to be labelled in the legend. This prevents
                the legend from containing labels for orbitals that have very
                little contribution in the plotting range.
            subplot (:obj:`bool`, optional): Plot the density of states for
                each element on seperate subplots. Defaults to ``False``.

        Returns:
            dict: The plotting data. Formatted with the following keys:

                "energies" (:obj:`numpy.ndarray`)
                    The energies.

                "mask" (:obj:`numpy.ndarray`)
                    A mask used to trim the density of states data and
                    prevent unwanted data being included in the output file.

                "lines" (:obj:`list`)
                    A :obj:`list` of :obj:`dict` containing the density data
                    and some metadata. Each line :obj:`dict` contains the keys:

                        "label" (:obj:`str`)
                            The label for the legend.

                        "dens" (:obj:`numpy.ndarray`)
                            The density of states data.

                        "colour" (:obj:`str`)
                            The colour of the line.

                        "alpha" (:obj:`float`)
                            The alpha value for line fill.

                "ymin" (:obj:`float`)
                    The minimum y-axis limit.

                "ymax" (:obj:`float`)
                    The maximum y-axis limit.
        """
        # mask needed to prevent unwanted data in pdf and for finding y limit
        dos = self._dos
        pdos = self._pdos
        mask = (dos.energies >= xmin - 0.05) & (dos.energies <= xmax + 0.05)
        plot_data = {'mask': mask, 'energies': dos.energies}
        spins = dos.densities.keys()
        ymax = 0

        if plot_total:
            lines = []
            tdos = {'label': 'Total DOS', 'dens': dos.densities,
                    'colour': 'k', 'alpha': 0.15}

            # subplot data formatted as a list of lists of dicts, with each
            # list of dicts being plotted on a seperate graph, if only one list
            # then solo plot
            lines.append([tdos])
            dmax = max([max(d[mask]) for d in dos.densities.values()])
            ymax = dmax if dmax > ymax else ymax
        elif not subplot:
            lines = [[]]  # need a blank list to add lines into

        # TODO: Fix broken behaviour if plot_total is off
        cutoff = (legend_cutoff / 100.) * (ymax / 1.05)

        for el, el_pdos in pdos.items():
            el_lines = []
            for orb in sort_orbitals(el_pdos):
                dmax = max([max(d[mask])
                            for d in el_pdos[orb].densities.values()])
                ymax = dmax if dmax > ymax else ymax
                label = None if dmax < cutoff else '{} ({})'.format(el, orb)
                colour = get_colour_for_element_and_orbital(el, orb, colours)
                el_lines.append({'label': label, 'alpha': 0.25,
                                 'colour': colour,
                                 'dens': el_pdos[orb].densities})
            if subplot:
                lines.append(el_lines)
            else:
                lines[0].extend(el_lines)

        ymax = ymax * empty_space / yscale
        ymin = 0 if len(spins) == 1 else -ymax
        plot_data.update({'lines': lines, 'ymax': ymax, 'ymin': ymin})
        return plot_data

    def get_plot(self, subplot=False, width=6., height=8., xmin=-6., xmax=6.,
                 yscale=1, colours=None, plot_total=True, legend_on=True,
                 num_columns=2, legend_frame_on=False, legend_cutoff=3,
                 xlabel='Energy (eV)', ylabel='Arb. units',
                 dpi=400, fonts=None, plt=None):
        """Get a :obj:`matplotlib.pyplot` object of the density of states.

        Args:
            subplot (:obj:`bool`, optional): Plot the density of states for
                each element on seperate subplots. Defaults to ``False``.
            width (:obj:`float`, optional): The width of the plot.
            height (:obj:`float`, optional): The height of the plot.
            xmin (:obj:`float`, optional): The minimum energy on the x-axis.
            xmax (:obj:`float`, optional): The maximum energy on the x-axis.
            yscale (:obj:`float`, optional): Scaling factor for the y-axis.
            colours (:obj:`dict`, optional): Use custom colours for specific
                element and orbital combinations. Specified as a :obj:`dict` of
                :obj:`dict` of the colours. For example::

                    {
                        'Sn': {'s': 'r', 'p': 'b'},
                        'O': {'s': '#000000'}
                    }

                The colour can be a hex code, series of rgb value, or any other
                format supported by matplotlib.
            plot_total (:obj:`bool`, optional): Plot the total density of
                states. Defaults to ``True``.
            legend_on (:obj:`bool`, optional): Plot the graph legend. Defaults
                to ``True``.
            num_columns (:obj:`int`, optional): The number of columns in the
                legend.
            legend_frame_on (:obj:`bool`, optional): Plot a frame around the
                graph legend. Defaults to ``False``.
            legend_cutoff (:obj:`float`, optional): The cut-off (in % of the
                maximum density of states within the plotting range) for an
                elemental orbital to be labelled in the legend. This prevents
                the legend from containing labels for orbitals that have very
                little contribution in the plotting range.
            xlabel (:obj:`str`, optional): Label/units for x-axis (i.e. energy)
            ylabel (:obj:`str`, optional): Label/units for y-axis (i.e. DOS)
            dpi (:obj:`int`, optional): The dots-per-inch (pixel density) for
                the image.
            fonts (:obj:`list`, optional): Fonts to use in the plot. Can be a
                a single font, specified as a :obj:`str`, or several fonts,
                specified as a :obj:`list` of :obj:`str`.
            plt (:obj:`matplotlib.pyplot`, optional): A
                :obj:`matplotlib.pyplot` object to use for plotting.

        Returns:
            :obj:`matplotlib.pyplot`: The density of states plot.
        """
        plot_data = self.dos_plot_data(yscale=yscale, xmin=xmin, xmax=xmax,
                                       colours=colours, plot_total=plot_total,
                                       legend_cutoff=legend_cutoff,
                                       subplot=subplot)

        if subplot:
            nplots = len(plot_data['lines'])
            plt = pretty_subplot(nplots, 1, width=width, height=height,
                                 dpi=dpi, plt=plt, fonts=fonts)
        else:
            plt = pretty_plot(width=width, height=height, dpi=dpi, plt=plt,
                              fonts=fonts)

        mask = plot_data['mask']
        energies = plot_data['energies'][mask]
        fig = plt.gcf()
        lines = plot_data['lines']
        spins = [Spin.up] if len(lines[0][0]['dens']) == 1 else \
            [Spin.up, Spin.down]

        for i, line_set in enumerate(plot_data['lines']):
            if subplot:
                ax = fig.axes[i]
            else:
                ax = plt.gca()

            for line, spin in itertools.product(line_set, spins):
                if spin == Spin.up:
                    label = line['label']
                    densities = line['dens'][spin][mask]
                elif spin == Spin.down:
                    label = ""
                    densities = -line['dens'][spin][mask]
                ax.fill_between(energies, densities, lw=0,
                                facecolor=line['colour'],
                                alpha=line['alpha'])
                ax.plot(energies, densities, label=label,
                        color=line['colour'], lw=line_width)

            ax.set_ylim(plot_data['ymin'], plot_data['ymax'])
            ax.set_xlim(xmin, xmax)

            ax.tick_params(axis='x', which='both', top='off')
            ax.tick_params(axis='y', which='both', labelleft='off',
                           labelright='off', left='off', right='off')

            loc = 'upper right' if subplot else 'best'
            ncol = 1 if subplot else num_columns
            if legend_on:
                ax.legend(loc=loc, frameon=legend_frame_on, ncol=ncol,
                          handlelength=2, prop={'size': label_size - 3})

        # no add axis labels and sort out ticks
        if subplot:
            ax.set_xlabel(xlabel, fontsize=label_size)
            fig.subplots_adjust(hspace=0)
            plt.setp([a.get_xticklabels() for a in fig.axes[:-1]],
                     visible=False)
            fig.text(0.08, 0.5, ylabel, fontsize=label_size, ha='left',
                     va='center', rotation='vertical', transform=ax.transAxes)
        else:
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)

        return plt


def get_colour_for_element_and_orbital(element, orbital, colours=None):
    """Get a colour for a particular elemental and orbital combination.

    If the element is not specified in the colours dictionary, a random colour
    will be used based on the list of 22 colours of maximum contast:
    http://www.iscc.org/pdf/PC54_1724_001.pdf

    This is cached in sumo.plotting.colour_cache and re-used for subsequent
    plots within this Python instance. To reset the cache, set
    ``sumo.plotting.colour_cache = {}``.

    Args:
        element (str): The element.
        orbital (str): The orbital.
        colours (:obj:`dict`, optional): Use custom colours for specific
            element and orbital combinations. Specified as a :obj:`dict` of
            :obj:`dict` of the colours. For example::

                {
                    'Sn': {'s': 'r', 'p': 'b'},
                    'O': {'s': '#000000'}
                }

            The colour can be a hex code, series of rgb value, or any other
            format supported by matplotlib.

    Returns:
        str: The colour.
    """

    def _get_colour_with_cache(element, orbital):
        """Return cached colour if available, or fetch and cache from cycle"""
        if element in colour_cache and orbital in colour_cache[element]:
            return colour_cache[element][orbital]
        else:
            colour = next(col_cycle)
            if element not in colour_cache:
                colour_cache[element] = {}
            colour_cache[element].update({orbital: colour})
            return colour

    if isinstance(colours, configparser.ConfigParser):
        try:
            return colours.get(element, orbital)
        except(configparser.NoSectionError, configparser.NoOptionError):
            return _get_colour_with_cache(element, orbital)

    elif isinstance(colours, dict):
        try:
            return colours[element][orbital]
        except(KeyError):
            return _get_colour_with_cache(element, orbital)

    elif colours is None:
        return _get_colour_with_cache(element, orbital)

    else:
        raise TypeError('Argument "colours" should be dict, '
                        'ConfigParser or None.')
