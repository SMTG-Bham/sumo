# coding: utf-8
# Copyright (c) Scanlon Materials Theory Group
# Distributed under the terms of the MIT License.

import logging

import numpy as np
import itertools as it

from matplotlib.ticker import MaxNLocator, AutoMinorLocator

from vaspy.plotting import pretty_plot, pretty_subplot, rgbline
from vaspy.electronic_structure.bandstructure import \
        get_projections_by_branches

from pymatgen.electronic_structure.plotter import BSPlotter
from pymatgen.electronic_structure.core import Spin

line_width = 1.5
empty_space = 1.05
label_size = 22
band_linewidth = 2


class VBSPlotter(BSPlotter):

    def __init__(self, bs):
        """Vaspy class for plotting band structures.

        This class is similar to the pymatgen BSPlotter class but overrides
        some methods to generate prettier plots.

        Further functionality, such as projected band structure plots are
        available.

        Args:
            bs (BandStructure): A pymatgen BandStructure object.
        """
        BSPlotter.__init__(self, bs)

    def get_plot(self, zero_to_efermi=True, ymin=-6., ymax=6.,
                 width=6., height=6., vbm_cbm_marker=False, dpi=400, plt=None,
                 dos_plotter=None, dos_options=None, dos_aspect=3, fonts=None):
        """Get a matplotlib object for the bandstructure plot.

        If spin polarised or metallic, blue lines are spin up, red lines are
        spin down. Otherwise, for semiconductors, blue lines indicate valence
        bands and orange indicates conduction bands.

        Args:
            zero_to_efermi (bool): Automatically subtract off the Fermi energy
                from the eigenvalues and plot (E-Ef).
            ymin (float): The y-axis (energy) minimum limit.
            ymax (float): The y-axis (energy) maximum limit.
            width (float): The width of the figure.
            height (float): The height of the figure.
            vbm_cbm_marker (bool): Plot markers to indicate the VBM and CBM
                locations.
            dpi (int): The dots-per-inch (pixel density) for the image.
            plt (pyplot object): Matplotlib pyplot object to use for plotting.
            fonts (list): List of fonts to use in the plot.
        """
        if dos_plotter:
            width = width + height/dos_aspect
            plt = pretty_subplot(1, 2, width, height, sharex=False, dpi=dpi,
                                 plt=plt, fonts=fonts,
                                 gridspec_kw={'width_ratios': [dos_aspect, 1],
                                              'wspace': 0})
            ax = plt.gcf().axes[0]
        else:
            plt = pretty_plot(width, height, dpi=dpi, plt=plt, fonts=fonts)
            ax = plt.gca()

        data = self.bs_plot_data(zero_to_efermi)
        dists = data['distances']
        eners = data['energy']

        if self._bs.is_spin_polarized or self._bs.is_metal():
            is_vb = True
        else:
            is_vb = self._bs.bands[Spin.up] <= self._bs.get_vbm()['energy']

        # nd is branch index, nb is band index, nk is kpoint index
        for nd, nb in it.product(range(len(data['distances'])),
                                 range(self._nb_bands)):
            e = eners[nd][str(Spin.up)][nb]

            # this check is very slow but works for now
            # colour valence bands blue and conduction bands orange
            if (self._bs.is_spin_polarized or self._bs.is_metal() or
                    np.all(is_vb[nb])):
                c = '#3953A4'
            else:
                c = '#FAA316'

            # plot band data
            ax.plot(dists[nd], e, ls='-', c=c, linewidth=band_linewidth)
            if self._bs.is_spin_polarized:
                e = eners[nd][str(Spin.down)][nb]
                ax.plot(dists[nd], e, 'r--', linewidth=band_linewidth)

        self._maketicks(ax)
        self._makeplot(ax, plt.gcf(), data, zero_to_efermi=zero_to_efermi,
                       vbm_cbm_marker=vbm_cbm_marker, width=width,
                       height=height, ymin=ymin, ymax=ymax,
                       dos_plotter=dos_plotter, dos_options=dos_options)
        plt.tight_layout()
        return plt

    def get_projected_rgb_plot(self, selection, zero_to_efermi=True,
                               ymin=-6., ymax=6., width=6., height=6.,
                               vbm_cbm_marker=False, dpi=400, plt=None,
                               dos_plotter=None, dos_options=None,
                               dos_aspect=3):
        """Get a matplotlib object for a projected rgb bandstructure plot.

        The band structure line color depends on the character of the band
        (either different elemental or orbital contributions). Each element/
        orbital is associated with red, green or blue and the corresponding rgb
        color depending on the character of the band is used. The method can
        only deal with up to 3 elements/orbitals.

        Spin up and spin down are differientiated by a '-' and a '--' line

        Args:
            selection: A list of tuples/strings identifying which elements
                and orbitals to project on to the band structure. These can be
                specified by both element and orbital, for example:

                    [('Bi', 's'), ('Bi', 'p'), ('S', 'p')]

                If just the element is specified then all the orbitals of
                that element are combined. For example:

                    [('Bi', 's'), ('Bi', 'p'), 'S']

                You can also choose to sum certain orbitals, by supplying a
                tuple of orbitals. For example:

                    [('Bi', 's'), ('Bi', 'p'), ('S', ('s', 'p', 'd'))]

                A maximum of 3 tuples can be plotted simultaneously (one for
                red, green and blue). The order of the tuples will affect which
                colour is used.
            zero_to_efermi (bool): Automatically subtract off the Fermi energy
                from the eigenvalues and plot (E-Ef).
            ymin (float): The y-axis (energy) minimum limit.
            ymax (float): The y-axis (energy) maximum limit.
            width (float): The width of the figure.
            height (float): The height of the figure.
            vbm_cbm_marker (bool): Plot markers to indicate the VBM and CBM
                locations.
            dpi (int): The dots-per-inch (pixel density) for the image.
            plt (pyplot object): Matplotlib pyplot object to use for plotting.
        """
        if dos_plotter:
            width = width + height/dos_aspect
            plt = pretty_subplot(1, 2, width, height, sharex=False, dpi=dpi,
                                 plt=plt,
                                 gridspec_kw={'width_ratios': [dos_aspect, 1],
                                              'wspace': 0})
            ax = plt.gcf().axes[0]
        else:
            plt = pretty_plot(width, height, dpi=dpi, plt=plt)
            ax = plt.gca()

        if len(selection) > 3:
            raise ValueError('Too many elements/orbitals specified (max 3)')

        data = self.bs_plot_data(zero_to_efermi)
        dists = data['distances']
        eners = data['energy']
        nbands = self._nb_bands
        nbranches = len(data['distances'])
        spins = self._bs.bands.keys()

        proj = get_projections_by_branches(self._bs, selection,
                                           normalise='select')

        # nd is branch index, nb is band index
        for spin, nd, nb in it.product(spins, range(nbranches), range(nbands)):
            colour = [proj[nd][i][spin][nb] for i in range(len(selection))]

            # if only two orbitals then just use red and blue
            if len(colour) == 2:
                colour.insert(1, np.zeros((len(dists[nd]))))

            ls = '-' if spin == Spin.up else '--'

            lc = rgbline(dists[nd], eners[nd][str(spin)][nb], colour[0],
                         colour[1], colour[2], alpha=1, linestyles=ls)
            ax.add_collection(lc)

            # alternative is to use scatter plot
            # ax.scatter(dists[nd], eners[nd][ns][nb], s=10, edgecolors='none',
            #            c=zip(colour[0], colour[1], colour[2],
            #                  np.ones(len(dists[nd]))))

        # plot the legend
        colours = ['r', 'g', 'b'] if len(selection) == 3 else ['r', 'b']
        for c, spec in zip(colours, selection):
            if type(spec) == str:
                label = spec
            else:
                label = '{} ({})'.format(spec[0], " + ".join(spec[1]))
            ax.scatter([-10000], [-10000], c=c, s=50, label=label)
        ax.legend(bbox_to_anchor=(0.95, 1), loc=2, frameon=False,
                  prop={'size': label_size-2}, handletextpad=0.1)

        # finish and tidy plot
        self._maketicks(ax)
        self._makeplot(ax, plt.gcf(), data, zero_to_efermi=zero_to_efermi,
                       vbm_cbm_marker=vbm_cbm_marker, ymin=ymin, ymax=ymax,
                       height=height, width=width,
                       dos_plotter=dos_plotter, dos_options=dos_options)
        return plt

    def _makeplot(self, ax, fig, data, zero_to_efermi=True,
                  vbm_cbm_marker=False, ymin=-6, ymax=6, height=6, width=6,
                  dos_plotter=None, dos_options=None):
        # draw line at Fermi level if not zeroing to e-Fermi
        if not zero_to_efermi:
            ef = self._bs.efermi
            ax.axhline(ef, linewidth=2, color='k')

        # set x and y limits
        ax.set_xlim(0, data['distances'][-1][-1])
        if self._bs.is_metal() and not zero_to_efermi:
            ax.set_ylim(self._bs.efermi + ymin, self._bs.efermi + ymax)
        else:
            ax.set_ylim(ymin, ymax)

        if vbm_cbm_marker:
            for cbm in data['cbm']:
                ax.scatter(cbm[0], cbm[1], color='#D93B2B', marker='o', s=100)
            for vbm in data['vbm']:
                ax.scatter(vbm[0], vbm[1], color='#0DB14B', marker='o', s=100)

        if dos_plotter:
            ax = fig.axes[1]
            dos_options.update({'xmin': ymin, 'xmax': ymax})
            self._makedos(ax, dos_plotter, dos_options)
        else:
            # keep correct aspect ratio square
            x0, x1 = ax.get_xlim()
            y0, y1 = ax.get_ylim()
            ax.set_aspect((height/width) * ((x1-x0)/(y1-y0)))

    def _makedos(self, ax, dos_plotter, dos_options):
        """This is basically the same as the VDOSPlotter get_plot function."""
        plot_data = dos_plotter.dos_plot_data(**dos_options)
        mask = plot_data['mask']
        energies = plot_data['energies'][mask]
        lines = plot_data['lines']
        spins = [Spin.up] if len(lines[0][0]['dens']) == 1 else \
            [Spin.up, Spin.down]
        for i, line_set in enumerate(plot_data['lines']):
            for line, spin in it.product(line_set, spins):
                if spin == Spin.up:
                    label = line['label']
                    densities = line['dens'][spin][mask]
                elif spin == Spin.down:
                    label = ""
                    densities = -line['dens'][spin][mask]
                ax.fill_betweenx(energies, densities, 0, lw=0,
                                 facecolor=line['colour'],
                                 alpha=line['alpha'])
                ax.plot(densities, energies, label=label,
                        color=line['colour'], lw=line_width)

            # x and y axis reversed versus normal dos plotting
            ax.set_ylim(dos_options['xmin'], dos_options['xmax'])
            ax.set_xlim(plot_data['ymin'], plot_data['ymax'])

            ax.tick_params(axis='y', which='both', top='off')
            ax.tick_params(axis='x', which='both', labelbottom='off',
                           labeltop='off', bottom='off', top='off')

            ax.legend(loc=2, frameon=False, ncol=1,
                      prop={'size': label_size - 3},
                      bbox_to_anchor=(1., 1.))

    def _maketicks(self, ax):
        """Utility method to add tick marks to a band structure."""
        # set y-ticks
        ax.yaxis.set_major_locator(MaxNLocator(6))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))

        # set x-ticks; only plot the unique tick labels
        ticks = self.get_ticks()
        unique_d = []
        unique_l = []
        if ticks['distance']:
            temp_ticks = list(zip(ticks['distance'], ticks['label']))
            unique_d.append(temp_ticks[0][0])
            unique_l.append(temp_ticks[0][1])
            for i in range(1, len(temp_ticks)):
                if unique_l[-1] != temp_ticks[i][1]:
                    unique_d.append(temp_ticks[i][0])
                    unique_l.append(temp_ticks[i][1])

        logging.info('Label positions:')
        for dist, label in list(zip(unique_d, unique_l)):
            logging.info('\t{:.4f}: {}'.format(dist, label))

        ax.set_xticks(unique_d)
        ax.set_xticklabels(unique_l)
        ax.xaxis.grid(True, c='k', ls='-', lw=line_width)
        ax.set_ylabel('Energy (eV)')
