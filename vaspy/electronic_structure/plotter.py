import itertools
import logging
import warnings
import vaspy.misc.plotting
import matplotlib
from matplotlib.ticker import MultipleLocator

from pymatgen.electronic_structure.plotter import BSPlotter
from pymatgen.electronic_structure.core import Spin


class VBSPlotter(BSPlotter):

    def __init__(self, bs):
        BSPlotter.__init__(self, bs)

    def _maketicks(self, plt):
        """Utility method to add tick marks to a band structure."""
        # set y-ticks
        ax = plt.gca()
        majorLocator = MultipleLocator(2)
        minorLocator = MultipleLocator(1)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_minor_locator(minorLocator)
        plt.tick_params(axis='both', which='major', labelsize=19)
        plt.tick_params(axis='both', which='minor', labelsize=19)

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
        ax.xaxis.grid(True, c='k', ls='-', lw=vaspy.misc.plotting._linewidth)

    def get_plot(self, zero_to_efermi=True, ymin=-6, ymax=6,
                 width=5, height=5, vbm_cbm_marker=False, plt=None,
                 dos_data=None):
        """
        Get a matplotlib object for the bandstructure plot.
        Blue lines are up spin, red lines are down
        spin.

        Args:
            zero_to_efermi: Automatically subtract off the Fermi energy from
                the eigenvalues and plot (E-Ef).
            ylim: Specify the y-axis (energy) limits; by default None let
                the code choose. It is vbm-4 and cbm+4 if insulator
                efermi-10 and efermi+10 if metal
        """
        if not plt:
            plt = vaspy.misc.plotting.pretty_plot(width, height, plt=plt)

        band_linewidth = 2

        data = self.bs_plot_data(zero_to_efermi)
        dists = data['distances']
        eners = data['energy']
        nkpts = len(dists[0])
        # nd is branch index, nb is branch index, nk is kpoint index
        for nd, nb in itertools.product(range(len(data['distances'])),
                                        range(self._nb_bands)):
            # colour valence bands blue and conduction bands orange in semiconds
            if (self._bs.is_spin_polarized or self._bs.is_metal() or
                    nb <= self._bs.get_vbm()['band_index'][Spin.up][0] + 1):
                c = '#3953A4'
            else:
                c = '#FAA316'

            e = [eners[nd][str(Spin.up)][nb][nk] for nk in range(nkpts)]
            plt.plot(dists[nd], e, ls='-', c=c, linewidth=band_linewidth)
            if self._bs.is_spin_polarized:
                e = [eners[nd][str(Spin.down)][nb][nk] for nk in range(nkpts)]
                plt.plot(dists[nd], e, 'r--', linewidth=band_linewidth)

        self._maketicks(plt)

        plt.ylabel('Energy (eV)')

        # draw line at Fermi level if not zeroing to e-Fermi
        if not zero_to_efermi:
            ef = self._bs.efermi
            plt.axhline(ef, linewidth=2, color='k')

        # set x and y limits
        plt.xlim(0, data['distances'][-1][-1])
        if self._bs.is_metal() and not zero_to_efermi:
            plt.ylim(self._bs.efermi + ymin, self._bs.efermi + ymax)
        else:
            plt.ylim(ymin, ymax)

        if vbm_cbm_marker:
            for cbm in data['cbm']:
                plt.scatter(cbm[0], cbm[1], color='#D93B2B', marker='o', s=100)
            for vbm in data['vbm']:
                plt.scatter(vbm[0], vbm[1], color='#0DB14B', marker='o', s=100)

        # keep correct aspect ratio square
        ax = plt.gca()
        x0, x1 = ax.get_xlim()
        y0, y1 = ax.get_ylim()
        ax.set_aspect((height/width) * ((x1-x0)/(y1-y0)))

        plt.tight_layout()
        return plt
