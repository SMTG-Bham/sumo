# coding: utf-8
# Copyright (c) Scanlon Materials Theory Group
# Distributed under the terms of the MIT License.

import logging

import numpy as np
import itertools as it

from scipy.interpolate import interp1d
from matplotlib.ticker import MaxNLocator, AutoMinorLocator

from vaspy.plotting import (pretty_plot, pretty_subplot, rgbline,
                            default_colours)
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

    def get_projected_plot(self, selection, mode='rgb', interpolate_factor=4,
                           circle_size=150, projection_cutoff=0.001,
                           zero_to_efermi=True, ymin=-6., ymax=6., width=6.,
                           height=6., vbm_cbm_marker=False, dpi=400, plt=None,
                           dos_plotter=None, dos_options=None,
                           dos_aspect=3):
        """Get a matplotlib object for a projected bandstructure plot.

        For the mode='rgb', spin up and spin down are differientiated by a '-'
        and a '--' line, otherwise spin up and spin down are plotted
        seperately.

        Args:
            selection (list): A list of tuples/strings identifying which
                elements and orbitals to project on to the band structure.
                These can be specified by both element and orbital, for example

                    [('Bi', 's'), ('Bi', 'p'), ('S', 'p')]

                If just the element is specified then all the orbitals of
                that element are combined. For example:

                    [('Bi', 's'), ('Bi', 'p'), 'S']

                You can also choose to sum certain orbitals, by supplying a
                tuple of orbitals. For example:

                    [('Bi', 's'), ('Bi', 'p'), ('S', ('s', 'p', 'd'))]

                If the plotting mode is 'rgb', a maximum of 3 orbital/element
                combinations can be plotted simultaneously (one for red, green
                and blue), otherwise an unlimited number of elements/orbitals
                can be selected.
            mode (str): Type of projected band structure to plot. Options are:
                "rgb": The band structure line color depends on the character
                    of the band. Each element/orbital contributes either red,
                    green or blue with the corresponding line colour a mixture
                    of all three colours. This mode only supports up to with up
                    to 3 elements/orbitals. The order of the tuples determines
                    which colour is used.
                "stacked": The element/orbital contributions are drawn as a
                    series of stacked circles, with the colour depending on
                    the composition of the band. The size of the circles can
                    be scaled using the stacked_marker_size option.
            circle_size (float): The size of the circles used in the 'stacked'
                plotting mode.
            projection_cutoff (float): Don't plot projections with normalised
                intensitites below this number. This option is useful for
                stacked plots, where small projections look nasty.
            interpolate_factor (int): The factor by which to interpolate the
                band structure (neccessary to make smooth lines). A larger
                number indicates greater interpolation.
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
        if mode == 'rgb' and len(selection) > 3:
            raise ValueError('Too many elements/orbitals specified (max 3)')
        elif mode == 'solo' and dos_plotter:
            raise ValueError('Solo mode plotting with DOS not supported')

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

        data = self.bs_plot_data(zero_to_efermi)
        nbranches = len(data['distances'])

        # Ensure we do spin up first, then spin down
        spins = sorted(self._bs.bands.keys(), key=lambda spin: -spin.value)

        proj = get_projections_by_branches(self._bs, selection,
                                           normalise='select')

        # nd is branch index
        for spin, nd in it.product(spins, range(nbranches)):

            # mask data to reduce plotting load
            bands = np.array(data['energy'][nd][str(spin)])
            mask = np.where(np.any(bands > ymin - 0.05, axis=1) &
                            np.any(bands < ymax + 0.05, axis=1))
            distances = data['distances'][nd]
            bands = bands[mask]
            weights = [proj[nd][i][spin][mask] for i in range(len(selection))]

            # interpolate band structure to improve smoothness
            dx = (distances[1] - distances[0]) / interpolate_factor
            temp_dists = np.arange(distances[0], distances[-1], dx)
            bands = interp1d(distances, bands, axis=1)(temp_dists)
            weights = interp1d(distances, weights, axis=2)(temp_dists)
            distances = temp_dists

            if mode == 'rgb':

                # colours aren't used now but needed later for legend
                colours = ['r', 'g', 'b']

                # if only two orbitals then just use red and blue
                if len(weights) == 2:
                    weights = np.insert(weights, 1, np.zeros(weights[0].shape),
                                        axis=0)
                    colours = ['r', 'b']

                ls = '-' if spin == Spin.up else '--'
                lc = rgbline(distances, bands, weights[0],
                             weights[1], weights[2], alpha=1, linestyles=ls)
                ax.add_collection(lc)

            elif mode == 'stacked':
                # TODO: Handle spin

                # use some nice custom colours first, then default colours
                colours = ['#3952A3', '#FAA41A', '#67BC47', '#6ECCDD',
                           '#ED2025']
                colours.extend(np.array(default_colours)/255)

                # very small cicles look crap
                weights[weights < projection_cutoff] = 0

                distances = list(distances) * len(bands)
                bands = bands.flatten()
                zorders = range(-len(weights), 0)
                for w, c, z in zip(weights, colours, zorders):
                    ax.scatter(distances, bands, c=c, s=circle_size*w**2,
                               zorder=z, rasterized=True)

        # plot the legend
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
