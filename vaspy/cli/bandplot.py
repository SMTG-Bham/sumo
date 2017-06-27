#!/usr/bin/env python
import os
import sys
import glob
import logging
import argparse

from vaspy.band import Symmetry, SymException, PathGen
from vaspy.vasp_input import Poscar
from vaspy.xmgrace import XmgraceDocument, XmgraceGraph
from vaspy.utils.bandgen import Bandgen, BandgenOptions


class BandPlotOptions(object):

    def __init__(self, options={}):
        if options:
            self.options = options
        else:
            parser = argparse.ArgumentParser(
                description="plot/(and generate) band structure and dos")
            parser.add_argument("prefix", help="prefix for generated files")
            parser.add_argument("-contcar", default="CONTCAR",
                                help="CONTCAR from which to determine "
                                "band path labels")

            symg = parser.add_mutually_exclusive_group()
            symg.add_argument("-tolerance", type=int, default=2,
                              choices=[1, 2, 3, 4, 5],
                              help="number of decimal places for symmetry "
                              "analysis tolerance - default: 2")
            symg.add_argument("-symmetry",
                              help="override kgen symmetry determination."
                              "should be in the form of 'centering,type', "
                              "e.g. 'P,orthorhombic'")

            parser.add_argument("-density", type=int, default=400,
                                help="kpoint line density - default: 400")
            parser.add_argument("-graphics", type=int,  default=0,
                                choices=[1, 2, 3],
                                help="generate graphics - 1 for pdf, 2 for png,"
                                " 3 for both")
            parser.add_argument("-vbm",
                                help="specify the valence band maximum if data "
                                "not available. write in the form of band_num,"
                                "energy - e.g -v 145,-34.62")
            parser.add_argument("-ignore", action="store_true",
                                help="ignore .dat files in current directory "
                                "and look for PROCARS (useful if you make a "
                                "mistake the first time)")
            parser.add_argument("-old", action="store_true",
                                help="Don't use the new path for orthorhombic P")
            parser.add_argument("-emin", type=float, default=-6,
                                help="Min energy in plot (relative to E_fermi)")
            parser.add_argument("-emax", type=float, default=6,
                                help="Max energy in plot (relative to E_fermi)")

            args = parser.parse_args()

            prefix = args.prefix+"_" if args.prefix else ""

            if args.vbm:
                band_num, energy = map(float, args.vbm.split(","))
                vbm = {'band': band_num, 'energy': energy}
            else:
                vbm = None

            self.options = {'prefix': prefix, 'tolerance': args.tolerance,
                            'symmetry': args.symmetry, 'density': args.density,
                            'graphics': args.graphics, 'contcar': args.contcar,
                            'vbm': vbm, 'ignore': args.ignore, 'old': args.old,
                            'emin': args.emin, 'emax': args.emax}


class BandPlot(object):

    def __init__(self, options):
        self.options = options.options
        self.band_data = self._load_bandstructures()
        self.axis_labels = self._get_band_labels()

    def _load_bandstructures(self):
        opts = {'prefix': self.options['prefix'], 'suffix': "PBE|HSE",
                'combine': False, 'search': False}
        bg_opts = BandgenOptions(opts)

        filenames = glob.glob('./*band*.dat')
        if filenames and not self.options['ignore']:
            filenames = [filename.strip('./') for filename in filenames]
            logging.info("found %d bandstructures in the current directory",
                         len(filenames))

            if not self.options['vbm']:
                vbm = {'energy': 0.0, 'band': float("-inf")}
                cbm = {'energy': 0.0, 'band': float("inf")}
                logging.info("\nno VBM or CBM data available "
                             "so will not colour bands")

        elif os.path.isfile("PROCAR"):
            bg = Bandgen(bg_opts)
            filenames = bg.generate_bandstructure()
            vbm = bg.vbm
            cbm = bg.cbm

        elif os.path.isfile("PROCAR1") or os.path.isfile("PROCAR01"):
            bg_opts.options['combine'] = True
            bg = Bandgen(bg_opts)
            filenames = bg.generate_bandstructure()
            vbm = bg.vbm
            cbm = bg.cbm

        else:
            bg_opts.options['combine'] = False
            bg_opts.options['search'] = True
            bg = Bandgen(bg_opts)
            filenames = bg.generate_bandstructure()
            vbm = bg.vbm
            cbm = bg.cbm

        filenames.sort()

        band_structures = []
        nkpoints = []
        for filename in filenames:
            band_path = []
            with open(filename) as f:
                lines = f.readlines()

            bands = []
            for line in lines:
                splits = line.split()
                if splits:
                    bands.append(map(float, splits))
                else:
                    band_path.append(bands)
                    bands = []

            nkpoints.append(len(band_path[0]))
            band_structures.append(band_path)

        nbands = len(band_structures[0])

        if self.options['vbm']:
            vbm = self.options['vbm']
            cbm = {'energy': 0.0, 'band': vbm['band'] + 1}

        return {'band_structures': band_structures, 'vbm': vbm, 'cbm': cbm,
                'nbands': nbands, 'nkpoints': nkpoints}

    def _get_band_labels(self):
        try:
            poscar = Poscar.load(self.options['contcar'])
        except:
            logging.info("unable to find %s, skipping label generation",
                         self.options['contcar'])
            return [{'labels': [], 'offsets': []}
                    for band in self.band_data['band_structures']]

        try:
            if self.options['symmetry'] is None:
                sym = Symmetry.from_poscar(poscar,
                                           tolerance=self.options['tolerance'])
                logging.info("\ndetermined symmetry to be: %s", sym)
            else:
                split_sym = self.options['symmetry'].split(",")
                if len(split_sym) != 2:
                    sys.exit("symmetry format invalid."
                             "Should be e.g. \"P,orthorhombic\"")

                centre, system = split_sym
                logging.info("specified symmetry is %s centred %s lattice",
                             centre, system)
                conv_cell = poscar.lattice_abc + poscar.lattice_angles
                sym = Symmetry({'spacegroup_symbol': centre,
                                'crystal_system': system}, conv_cell)

            paths = PathGen(poscar, sym, density=self.options['density'],
                    old=self.options['old'])

        except SymException as e:
            print e.message
            logging.info("skipping label generation")
            return [{'labels': [], 'offsets': []}
                    for band in self.band_data['band_structures']]

        try:
            label_data = []
            for path, step_size, nkpoints in zip(paths.labels,
                                                 paths.nkpoints_per_step,
                                                 self.band_data['nkpoints']):
                logging.info("for graph %d the labels are: %s",
                             paths.labels.index(path)+1, " -> ".join(path))

                calc_kpoints = sum(step_size) + 1
                if nkpoints != calc_kpoints:
                    raise
                label_data.append({'labels': path, 'offsets': step_size})

        except Exception as e:
            logging.info("\ngenerated kpoint path does not have same number of "
                         "kpoints (%d) as the band structure supplied (%d)"
                         "therefore can't apply labels",
                         calc_kpoints, nkpoints)
            return [{'labels': [], 'offsets': []}
                    for band in self.band_data['band_structures']]

        return label_data

    def to_agr(self):
        doc = XmgraceDocument(True)

        for struct, label_data, kpoints in\
                zip(self.band_data['band_structures'], self.axis_labels,
                    self.band_data['nkpoints']):

            struct_id = self.band_data["band_structures"].index(struct)

            # Use this to know where to put y-axis labels
            if struct_id == 0:
                pos = 0
            elif struct_id == len(self.band_data["band_structures"]) - 1:
                pos = 1
            else:
                pos = 2

            graph = XmgraceGraph(pos=pos, num_kpoints=kpoints,
                                 emin=self.options['emin'],
                                 emax=self.options['emax'])

            for band in struct:
                normalised = [[kpoint, energy - self.band_data['vbm']['energy']]
                              for kpoint, energy in band]
                if struct.index(band) <= (self.band_data['cbm']['band'] - 2):
                    color = 4
                else:
                    color = 11
                graph.add_xy(normalised, color, 2)

            axis_data = []
            total = 1
            if label_data['labels']:
                # hack to make the lists the same length
                label_data['offsets'].append(0)
                for label, offset in zip(label_data['labels'],
                                         label_data['offsets']):
                    label = '\\xG' if label == "GP" else label
                    axis_data.append([label, total])
                    total += offset

            graph.set_axis_labels(axis_data)
            doc.add_graph(graph)

        filename = self.options['prefix']+"band.agr"
        doc.to_file(filename)

        return filename

    @staticmethod
    def cmd_exists(cmd):
        return any(
            os.access(os.path.join(path, cmd), os.X_OK)
            for path in os.environ["PATH"].split(os.pathsep)
        )

    def to_graphic(self, agr_file):
        graphics = self.options['graphics']

        if self.cmd_exists("xmgrace"):
            xmgrace = "xmgrace"
        elif self.cmd_exists("/shared/ucl/apps/xmgrace/grace/bin/xmgrace"):
            xmgrace = "/shared/ucl/apps/xmgrace/grace/bin/xmgrace"  # legion

        if not xmgrace:
            logging.error("\ncan't find xmgrace therefore can't make graphics")
            return

        logging.info("\ngenerating eps, xmgrace path is: %s", xmgrace)

        filename = self.options['prefix'] + "band.eps"
        command = "%s %s -nosafe -hardcopy -hdevice EPS -printfile %s" % \
            (xmgrace, agr_file, filename)

        os.system(command)

        if graphics == 1 or graphics == 3:
            pdf_file = self.options['prefix'] + "band.pdf"
            if self.cmd_exists("epstopdf"):
                logging.info("writing pdf file")
                os.system("epstopdf %s" % filename)
            elif self.cmd_exists("convert"):
                logging.info("epstopdf not available so pdf won't look great")
                os.system("convert -density 400 %s %s" % (filename, pdf_file))
            else:
                logging.info("no software available to convert eps to pdf")

        if graphics == 2 or graphics == 3:
            png_file = self.options['prefix'] + "band.png"
            if not self.cmd_exists("convert"):
                logging.info("no software available to convert eps to png")
            else:
                logging.info("writing png file")
                os.system("convert -density 400 %s %s" % (filename, png_file))

def main():
    logging.basicConfig(filename='bandplot.log', level=logging.DEBUG,
                        format='%(message)s', filemode='w')
    console = logging.StreamHandler()
    logging.getLogger('').addHandler(console)
    options = BandPlotOptions()
    bdplot = BandPlot(options)
    agr_file = bdplot.to_agr()
    bdplot.to_graphic(agr_file)
    sys.exit(0)

if __name__ == '__main__':
    main()
