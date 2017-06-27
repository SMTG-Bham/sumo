import sys
import re
import os
import mmap
import logging
import numpy as np

from xml.etree.cElementTree import iterparse
from collections import defaultdict
from operator import itemgetter
from itertools import izip

from vaspy.vasp_input import Poscar, Incar
from vaspy.electronic_structure import VaspDos, Dos


class Procar:

    def __init__(self, kpoints, num_bands, num_ions, vbm, cbm, hybrid):
        self._hybrid = hybrid
        self._kpoints = []
        # remove weighted kpoints if the calculation hybrid
        if hybrid:
            for kpoint in kpoints:
                if kpoint['weight'] == 0:
                    self._kpoints.append(kpoint)
        else:
            self._kpoints = kpoints
        self._num_bands = num_bands
        self._num_ions = num_ions
        self._vbm = vbm
        self._cbm = cbm
        self._eigenvals = self._calc_eigenvals()

    def _calc_eigenvals(self):
        eigenvals = []
        for band in range(self._num_bands):
            kpoint_array = []
            for kpoint in range(len(self._kpoints)):
                if self.hybrid and self._kpoints[kpoint]['weight'] > 0:
                    continue
                kpoint_array.append(self._kpoints[kpoint]['bands']
                                                 [band]['energy'])
            eigenvals.append(kpoint_array)
        return eigenvals

    def str_VBM(self):
        return "VBM has an energy of %.8f at k-point %s (#%d)" % \
            (self.vbm['energy'], self.vbm['kpoint_coords'], self.vbm['kpoint'])

    def str_CBM(self):
        return "CBM has an energy of %.8f at k-point %s (#%d)" % \
            (self.cbm['energy'], self.cbm['kpoint_coords'], self.cbm['kpoint'])

    def str_bandgap(self):
        return "the band gap energy is %.3f eV" % self.bandgap

    @staticmethod
    def load(filename="PROCAR"):
        try:
            procar_file = open(filename)
            with procar_file as f:
                return Procar.from_file(f)
        except IOError:
            sys.exit("Cannot find " + filename)
        except (ValueError, IndexError):
            try:
                # This is a horrific hack, but i cba fixing it
                procar_file = open(filename)
                with procar_file as f:
                    return Procar.from_file(f, gapy_procar=True)
            except ValueError:
                procar_file = open(filename)
                with procar_file as f:
                    return Procar.from_file(f, zhenyu=True)

    @staticmethod
    def from_file(f, gapy_procar=False, zhenyu=False):
        f.readline()  # skip header line
        matches = re.search('[^\d]*([\d]*)[^\d]*([\d]*)[^\d]*([\d]*)',
                            f.readline())
        num_kpoints, num_bands, num_ions = map(int, matches.group(1, 2, 3))

        f.readline()

        vbm = {'kpoint': 0, 'band': 0, 'energy': float('-inf')}
        cbm = {'kpoint': 0, 'band': 0, 'energy': float('inf')}

        skipped = 0
        kpoints = []
        soc = 0

        for kpoint in range(num_kpoints):
            kpoint_line = f.readline()
            # can't use split cos sometimes it comes out as, "0.00-0.33"
            coords = map(float, [kpoint_line[18:28],
                                 kpoint_line[29:39],
                                 kpoint_line[40:50]])
            weight = float(kpoint_line[65:74])

            if weight > 0:
                skipped += 1

            f.readline()

            bands = []
            band_line = f.readline()

            for i in range(num_bands):
                band_num = i + 1
                band_e, band_occ = map(float, itemgetter(4, 7)
                                       (band_line.split()))

                if band_e > vbm['energy'] and band_occ >= 0.9:
                    vbm = {'kpoint': kpoint+1, 'kpoint_coords': coords,
                           'energy': band_e, 'band': band_num}

                if (band_e < cbm['energy'] or vbm['band'] > cbm['band'])\
                        and band_num > vbm['band'] and band_occ <= 0.10:
                    cbm = {'kpoint': kpoint+1, 'kpoint_coords': coords,
                           'energy': band_e, 'band': band_num}

                band = {'energy': band_e, 'occupancy': band_occ}
                bands.append(band)

                # find out if spin orbit coupling is on, unfortunately can't use
                # the occupancy of the first band as spin polarised calcs also
                # have a occ of 1 so have to check with trial and error
                if soc == 0:
                    p = 0
                    for _ in range(num_ions+4):
                        f.readline()
                    band_line = f.readline()
                    if band_line[0] == 'b':
                        soc = 1  # No SOC
                        gaps = 0
                    else:
                        band_line = f.readline()
                        if band_line[0] == 'b':
                            soc = 1  # No SOC
                            gaps = 1  # the system has f orbitals
                        else:
                            gaps = 4 if gapy_procar else 0
                            soc = 4  # SOC is turned on so now we to move ahead
                            n = 3
                            if gapy_procar:
                                logging.info("Weird PROCAR detected: "
                                             "Hello James!")
                            if zhenyu:
                                logging.info("Weird PROCAR detected: "
                                             "Hello Zhenyu!")
                                n = 2
                                soc = 3
                                p = 1
                            logging.info("   spin orbit coupling detected")

                            for _ in range(((num_ions - p) * n) + 1 + gaps):
                                f.readline()
                            band_line = f.readline()
                else:
                    for _ in range((soc * num_ions) + p + (soc - p*2) + 3 + gaps):
                        f.readline()
                    band_line = f.readline()

            kpoints.append({'coords': coords, 'weight': weight, 'bands': bands})

        if skipped == num_kpoints:
            print "   non hybrid PROCAR detected, not skipping any k-points"
            hybrid = False
        else:
            print "   hybrid PROCAR detected, skipping %d kpoints" % skipped
            hybrid = True
            num_kpoints += -skipped  # why did i do this ;( ???

        return Procar(kpoints, num_bands, num_ions, vbm, cbm, hybrid)

    def to_file(self, filename):
        with open(filename, 'w') as f:
            for band in range(self.num_bands):
                for kpoint in range(self.num_kpoints):
                    line = "%d  %.8f\n" %\
                        (kpoint+1, self._eigenvals[band][kpoint])
                    f.write(line)
                f.write("\n")

    @staticmethod
    def find_procars(search=True):
        if search:
            logging.info("searching for matching directories")
            dirs = next(os.walk('.'))[1]
            r = re.compile('band_(\d*)_(HSE|PBE)', re.I)
            r_split = re.compile('band_(\d*)_split_(\d*)_(HSE|PBE)', re.I)
            nosplit_band_dirs = filter(r.match, dirs)
            split_band_dirs = filter(r_split.match, dirs)

            if nosplit_band_dirs and split_band_dirs:
                sys.exit("both split and non split band paths found, "
                         "too confusing")
            elif not nosplit_band_dirs and not split_band_dirs:
                logging.info("no matching folders found")

        if nosplit_band_dirs or split_band_dirs:
            band_dirs = split_band_dirs + nosplit_band_dirs
            band_dirs.sort()

            filenames = []
            for folder in band_dirs:
                filename = folder + "/PROCAR"
                if os.path.isfile(filename):
                    logging.info("found PROCAR in %s", folder)
                else:
                    sys.exit("no PROCAR in %s" % folder)

                # Organise the filenames so we can use them more easily later
                num = int(folder.split("_")[1])
                if len(filenames) >= num:
                    filenames[num-1].append(filename)
                else:
                    filenames.append([filename])

            logging.info("\n%d band structure paths found from a total of %d "
                         "PROCARs", len(filenames), sum(map(len, filenames)))
        else:
            if os.path.isfile("PROCAR"):
                filenames = [["PROCAR"]]
                logging.info("loaded PROCAR")
            else:
                sys.exit("no PROCAR in current directory")
        return [Procar.combine(files) for files in filenames]

    @property
    def eigenvals(self):
        return self._eigenvals

    @property
    def num_kpoints(self):
        return len(self.eigenvals[0])

    @property
    def kpoints(self):
        kpoints = []
        for kpoint in self._kpoints:
            kpoints.append(kpoint['coords'])
        return kpoints

    @property
    def num_bands(self):
        return self._num_bands

    @property
    def num_ions(self):
        return self._num_ions

    @property
    def cbm(self):
        return self._cbm

    @property
    def vbm(self):
        return self._vbm

    @property
    def hybrid(self):
        return self._hybrid

    @property
    def bandgap(self):
        return (self.cbm['energy'] - self.vbm['energy'])

    @staticmethod
    def combine(filenames):
        kpoints = []
        vbm = {'energy': float('-inf')}
        cbm = {'energy': float('inf')}
        for f in filenames:
            procar = Procar.load(f)
            kpoints += procar._kpoints
            if procar.vbm['energy'] > vbm['energy']:
                vbm = procar.vbm
            if procar.cbm['energy'] < cbm['energy']:
                cbm = procar.cbm
        num_bands = len(kpoints[0]['bands'])
        num_ions = procar.num_ions
        hybrid = procar.hybrid
        return Procar(kpoints, num_bands, num_ions, vbm, cbm, hybrid)


# this class is currently flithy and needs sorting out, but it works
class Outcar:

    def __init__(self, outcar_file):
        self.filename = outcar_file
        with open(outcar_file) as f:
            s = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
            try:
                fermi_matches = re.findall('E-fermi\s*:\s+(\d*\.\d*)\s*.*', s)
                self.efermi = float(fermi_matches[-1])
            except:
                logging.error("Can't find E-fermi")

            try:
                energy_matches = re.findall('.*entropy=.*=\s*(-*\d*\.\d*).*', s)
                self.total_energy = float(energy_matches[-1])
            except:
                logging.error("Can't find total energy")

    @staticmethod
    def load(filename="OUTCAR", die=True):
        # not sure what relies on this method at the moment so leaving
        # functionality the same
        try:
            open(filename)
        except IOError:
            print "Cannot find " + filename
            if die:
                sys.exit(1)
        return Outcar(filename)

    def read_core_state_eigen(self):
        """
        Read the core state eigenenergies at each ionic step.
        This has been taken from the Outcar class of pymatgen.

        Returns:
            A list of dict over the atom such as [{"AO":[core state eig]}].
            The core state eigenenergy list for each AO is over all ionic
            step.

        Example:
            The core state eigenenergy of the 2s AO of the 6th atom of the
            structure at the last ionic step is [5]["2s"][-1]

        """

        with open(self.filename, "rt") as foutcar:
            line = foutcar.readline()
            while line != "":
                line = foutcar.readline()
                if "NIONS =" in line:
                    natom = int(line.split("NIONS =")[1])
                    cl = [defaultdict(list) for i in range(natom)]
                if "the core state eigen" in line:
                    iat = -1
                    while line != "":
                        line = foutcar.readline()
                        # don't know number of lines to parse without knowing
                        # specific species, so stop parsing when we reach
                        # "E-fermi" instead
                        if "E-fermi" in line:
                            break
                        data = line.split()
                        # data will contain odd number of elements if it is
                        # the start of a new entry, or even number of elements
                        # if it continues the previous entry
                        if len(data) % 2 == 1:
                            iat += 1  # started parsing a new ion
                            data = data[1:]  # remove element with ion number
                        for i in range(0, len(data), 2):
                            cl[iat][data[i]].append(float(data[i + 1]))
        return cl

    def read_last_core_state_eigen(self):
        """
        Read the core state eigenenergies from the last ionic step.

        Returns:
            A list of dict over the atom such as [{"AO":[core state eig]}].
            The core state eigenenergy list for each AO is over all ionic
            step.

        Example:
            The core state eigenenergy of the 2s AO of the 6th atom of the
            structure is [5]["2s"]

        """

        with open(self.filename, "rt") as foutcar:
            line = foutcar.readline()
            while line != "":
                line = foutcar.readline()
                if "NIONS =" in line:
                    natom = int(line.split("NIONS =")[1])
                    cl = [defaultdict(list) for i in range(natom)]
                if "the core state eigen" in line:
                    iat = -1
                    while line != "":
                        line = foutcar.readline()
                        # don't know number of lines to parse without knowing
                        # specific species, so stop parsing when we reach
                        # "E-fermi" instead
                        if "E-fermi" in line:
                            break
                        data = line.split()
                        # data will contain odd number of elements if it is
                        # the start of a new entry, or even number of elements
                        # if it continues the previous entry
                        if len(data) % 2 == 1:
                            iat += 1  # started parsing a new ion
                            data = data[1:]  # remove element with ion number
                        for i in range(0, len(data), 2):
                            cl[iat][data[i]] = float(data[i + 1])
        return cl


class Doscar:
    """
    Class to read in and process a DOSCAR file.

    There is no Doscar object as instead the electronic_structure.VaspDos
    should be used.

    The way vasp prints the DOS data depends on which parameters were used
    for the calculation. If LSORBIT = 0/1/2/10, the data will appear as:

        E  s  p  d    (4 colummns)

    Or if an atom has f orbitals, these will appear on the next line:

        E  s  p  d    (4 colummns)
        f             (1 column)

    If NSPIN = 2 is set then the orbitals will be split into the up and down
    contributions resulting:

        E  s(up)  s(down)  p(up)  p(down)  d(up)  d(down)    (7 columns)

    Or 7 columns then 2 columns on a new line if an atom has f orbitals.

    If LORBIT = 11, the orbitals will be further split into their contributions:

        E  s  py  pz  px  dxy  dyz  dz2  dxz  dx2-y2    (10 columns)

    Or 10 columns then 7 columns on a new line if an atom has f orbitals.

    If SOC is turned on with LSORBIT = .TRUE., then the data will be provided
    as follows, with a total of 13 columns:

        E  s(total) s(mx)  s(my)  s(mz)  p(total)  p(mx)  p(my)  p(mz)  [...]

    Or 13 columns then 4 columns on a new line if an atom has f orbitals.

    When NSPIN = 2 and LORBIT = 11, the data will split into 19 columns:

        E  s(up)  s(down)  py(up)  py(down)  pz(up)  pz(down)  px(up)  [...]

    Or 17 columns then 14 columns on a new line if an atom has f orbitals.

    If LSORBIT is set and LORBIT = 11 then the data will be printed with a
    total of 37 columns as:

        E  s(total)  s(mx)  s(my)  s(mz)  py(total)  py(mx)  py(my)  [...]

    Or 37 columns then 28 columns on a new line if an atom has f orbitals.

    Knowing this makes it possible to tell from the DOSCAR file which parameters
    were set in the VASP run.
    """

    @staticmethod
    def load(filename, find_poscar=True, specific_atoms=None):
        try:
            path = os.path.dirname(filename)
            poscar = Poscar.load(os.path.join(path, "CONTCAR"))
            species_info = zip(poscar.atomic_symbols, poscar.natoms)
        except Exception as e:
            logging.info("unable to load species info from CONTCAR/POSCAR\n" +
                         "cannot continue")
            print e.message
            sys.exit()

        with open(filename, "rt") as f:
            return Doscar.from_file(f, species_info,
                                    specific_atoms=specific_atoms)

    @staticmethod
    def from_file(f, species_info, specific_atoms=None):
        for _ in range(5):
            line = f.readline()  # skip header lines
        if not line.strip():
            raise ValueError("empty DOSCAR, check LORBIT is set")

        info_line = f.readline().split()
        nedos = int(info_line[2])
        efermi = float(info_line[3])
        logging.info("e-fermi is: %.3f" % efermi)

        data = _parse_tarray(f, nedos)
        nspin = 1
        if len(data[0]) > 3:
            logging.info("spin polarized detected")
            nspin = 2

        tdos = {}
        idos = {}
        energies = data[:, 0]
        for spin in range(nspin):
            # for nspin = 1, col = 1, 2 for nspin = 2, col = 1, 3 -> 2, 4
            tdos[spin] = data[:, 1 + spin]
            idos[spin] = data[:, 1 + spin + nspin]

        f.readline()  # skip the info line
        lm = None
        soc = None
        f_orbs = None
        pdoss = {}
        for el, natoms in species_info:
            el_pdoss = []
            for n in range(natoms):
                # work out the type of dos that has been run
                if not f_orbs:
                    data = _parse_tarray(f, nedos)
                else:
                    data = _parse_tarray(f, nedos*2)
                if soc is None:
                    l = len(data[0])
                    f_orbs = True if l != len(data[1]) else False
                    soc = True if (l == 13 or l == 37) else False  # skip columns
                    lm = True if (l == 10 or l == 19 or l == 37) else False
                    if soc == True:
                        logging.info("spin orbit coupling detected")
                    if lm:
                        logging.info("decomposed doscar detected")
                    if f_orbs:
                        logging.info("f orbitals detected")
                        # have to do this to correct for our bad inital parsing
                        data = np.append(data, _parse_tarray(f, nedos))

                pdos = Doscar.read_pdos(data, lm, f_orbs, nspin, soc)
                el_pdoss.append(pdos)
                f.readline()  # skip info line
            pdoss[el] = el_pdoss
        logging.info("finished loading DOSCAR")
        return VaspDos(efermi, energies, tdos, pdoss, dict(species_info))

    @staticmethod
    def read_pdos(data, lm, f, nspin, soc):
        """
        Reads in a section of pdos data for a particular atom.

        Args:
            data: A numpy array of all the orbitals spat out according to vasp
            lm: True if the DOSCAR is lm-decomposed
            f: True if the atom has f orbitals
            orbs: The orbitals to populate with data
            nspin: 2 for spin polarised 1 if not
            soc: The offset used if SOC is turned on

        Returns:
            A dict containing the densities in the form {Orbital:{Spin:Dens}}
        """
        if f:
            # if there are f orbitals the DOSCAR looks like:
            # energy s p d
            # f
            # energy s p d ... etc
            # so we need to pair up each two lines
            a = iter(data)
            pairs = izip(a, a)
            data = np.array([np.append(a, b) for a, b in pairs])
        if lm:
            orbs = Dos.lm_f_orbitals if f else Dos.lm_orbitals
        else:
            orbs = Dos.f_orbitals if f else Dos.orbitals
        col_skip = 4 if soc else 1
        pdos = defaultdict(dict)
        col = 1
        for orb in orbs:
            for spin in range(nspin):
                pdos[orb][spin] = data[:, col]
                col += col_skip
        return pdos


class Vasprun():

    def __init__(self, filename):
        self.filename = filename
        with open(filename, "rt") as f:
            self._parse(f)
        self.nionic_steps = len(self.ionic_steps)

    def _parse(self, stream):
        self.efermi = None
        self.other_dielectric = {}
        ionic_steps = []
        parsed_header = False
        for event, elem in iterparse(stream):
            tag = elem.tag
            if not parsed_header:
                if tag == "generator":
                    self.generator = self._parse_params(elem)
                elif tag == "incar":
                    self.incar = self._parse_params(elem)
                elif tag == "parameters":
                    self.parameters = self._parse_params(elem)
            if tag == "calculation":
                parsed_header = True
                ionic_steps.append(self._parse_calculation(elem))
            elif tag == "dielectricfunction":
                if ("comment" not in elem.attrib) or \
                   elem.attrib["comment"] == "INVERSE MACROSCOPIC DIELECTRIC TENSOR (including local field effects in RPA (Hartree))":
                    self.dielectric = self._parse_diel(elem)
                else:
                    self.other_dielectric[elem.attrib["comment"]] = \
                        self._parse_diel(elem)
        self.ionic_steps = ionic_steps
        self.vasp_version = self.generator["version"]

    @property
    def converged(self):
        """
        Returns true if a relaxation run is converged.
        """
        return self.converged_electronic and self.converged_ionic

    @property
    def converged_electronic(self):
        """
        Checks that electronic step convergence has been reached in the final
        ionic step
        """
        final_esteps = self.ionic_steps[-1]["electronic_steps"]
        if 'LEPSILON' in self.incar and self.incar['LEPSILON']:
            i = 1
            to_check = set(['e_wo_entrp', 'e_fr_energy', 'e_0_energy'])
            while set(final_esteps[i].keys()) == to_check:
                i += 1
            return i + 1 != self.parameters["NELM"]
        return len(final_esteps) < self.parameters["NELM"]

    @property
    def converged_ionic(self):
        """
        Checks that ionic step convergence has been reached, i.e. that vasp
        exited before reaching the max ionic steps for a relaxation run
        """
        nsw = self.parameters.get("NSW", 0)
        return nsw <= 1 or len(self.ionic_steps) < nsw

    def _parse_params(self, elem):
        params = {}
        for c in elem:
            name = c.attrib.get("name")
            if c.tag not in ("i", "v"):
                p = self._parse_params(c)
                if name == "response functions":
                    # Delete duplicate fields from "response functions",
                    # which overrides the values in the root params.
                    p = {k: v for k, v in p.items() if k not in params}
                params.update(p)
            else:
                ptype = c.attrib.get("type")
                val = c.text.strip() if c.text else ""
                if c.tag == "i":
                    params[name] = _parse_parameters(ptype, val)
                else:
                    params[name] = _parse_v_parameters(ptype, val,
                                                       self.filename, name)
        elem.clear()
        return Incar(params)

    def _parse_calculation(self, elem):
        try:
            istep = {i.attrib["name"]: float(i.text)
                     for i in elem.find("energy").findall("i")}
        except AttributeError:  # not all calculations have an energy
            istep = {}
            pass
        esteps = []
        for scstep in elem.findall("scstep"):
            try:
                d = {i.attrib["name"]: _vasprun_float(i.text)
                     for i in scstep.find("energy").findall("i")}
                esteps.append(d)
            except AttributeError:  # not all calculations have an energy
                pass
        for va in elem.findall("varray"):
            istep[va.attrib["name"]] = _parse_varray(va)
        istep["electronic_steps"] = esteps
        elem.clear()
        return istep

    @property
    def final_energy(self):
        """
        Final energy from the vasp run.
        """
        try:
            return self.ionic_steps[-1]["electronic_steps"][-1]["e_0_energy"]
        except (IndexError, KeyError):
            # not all calculations have a total energy, i.e. GW
            return np.inf

#    def _parse_dos(self, elem):
#        efermi = float(elem.find("i").text)
#        energies = None
#        tdensities = {}
#        idensities = {}
#
#        for s in elem.find("total").find("array").find("set").findall("set"):
#            data = np.array(_parse_varray(s))
#            energies = data[:, 0]
#            spin = 0 if s.attrib["comment"] == "spin 1" else 1
#            tdensities[spin] = data[:, 1]
#            idensities[spin] = data[:, 2]
#
#        pdoss = []
#        partial = elem.find("partial")
#        if partial is not None:
#            orbs = [ss.text for ss in partial.find("array").findall("field")]
#            orbs.pop(0)
#            lm = any(["x" in s for s in orbs])
#            for s in partial.find("array").find("set").findall("set"):
#                pdos = {}
#
#                for ss in s.findall("set"):
#                    spin = 0 if ss.attrib["comment"] == "spin 1" else 1
#                    data = np.array(_parse_varray(ss))
#                    nrow, ncol = data.shape
#                    for j in range(1, ncol):
#                        if lm:
#                            orb = Doscar.orbital_parts[j - 1]
#                        else:
#                            orb = orbs[j - 1].strip().upper()
#                        pdos[orb][spin] = data[:, j]
#                pdoss.append(pdos)
#        elem.clear()
#        return VaspDos(energies, tdensities, idensities, pdoss)


def _parse_varray(elem):
    return [[float(i) for i in v.text.split()] for v in elem]


def _parse_tarray(f, nlines):
    data = []
    for i in range(nlines):
        data.append([float(i) for i in f.readline().split()])
    return np.array(data)


def _parse_parameters(val_type, val):
    """
    Helper function to convert a Vasprun parameter into the proper type.
    Boolean, int and float types are converted.
    Args:
        val_type: Value type parsed from vasprun.xml.
        val: Actual string value parsed for vasprun.xml.
    """
    if val_type == "logical":
        return val == "T"
    elif val_type == "int":
        return int(val)
    elif val_type == "string":
        return val.strip()
    else:
        return float(val)


def _parse_v_parameters(val_type, val, filename, param_name):
    """
    Helper function to convert a Vasprun array-type parameter into the proper
    type. Boolean, int and float types are converted.
    Args:
        val_type: Value type parsed from vasprun.xml.
        val: Actual string value parsed for vasprun.xml.
        filename: Fullpath of vasprun.xml. Used for robust error handling.
            E.g., if vasprun.xml contains \*\*\* for some Incar parameters,
            the code will try to read from an INCAR file present in the same
            directory.
        param_name: Name of parameter.
    Returns:
        Parsed value.
    """
    if val_type == "logical":
        val = [i == "T" for i in val.split()]
    elif val_type == "int":
        try:
            val = [int(i) for i in val.split()]
        except ValueError:
            # Fix for stupid error in vasprun sometimes which displays
            # LDAUL/J as 2****
            val = _parse_from_incar(filename, param_name)
            if val is None:
                raise IOError("Error in parsing vasprun.xml")
    elif val_type == "string":
        val = val.split()
    else:
        try:
            val = [float(i) for i in val.split()]
        except ValueError:
            # Fix for stupid error in vasprun sometimes which displays
            # MAGMOM as 2****
            val = _parse_from_incar(filename, param_name)
            if val is None:
                raise IOError("Error in parsing vasprun.xml")
    return val


def _vasprun_float(f):
    """
    Large numbers are often represented as ********* in the vasprun.
    This function parses these values as np.nan
    """
    try:
        return float(f)
    except ValueError as e:
        f = f.strip()
        if f == '*' * len(f):
            warnings.warn('Float overflow (*******) encountered in vasprun')
            return np.nan
        raise e
