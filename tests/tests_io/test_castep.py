import json
import os
import unittest

try:
    from importlib.resources import files as ilr_files
except ImportError:  # Python < 3.9
    from importlib_resources import files as ilr_files

from monty.io import gzip
from monty.json import MontyDecoder
from numpy.testing import assert_array_almost_equal
from pymatgen.core.structure import Structure
from pymatgen.electronic_structure.core import Spin

from sumo.io.castep import (
    CastepCell,
    CastepPhonon,
    labels_from_cell,
    read_bands_eigenvalues,
    read_bands_header,
    read_dos,
)

_ry_to_ev = 13.605693009


class CastepCellTestCase(unittest.TestCase):
    def setUp(self):
        self.si_cell = os.path.join(ilr_files("tests"), "data", "Si", "Si2.cell")
        self.si_cell_alt = os.path.join(
            ilr_files("tests"), "data", "Si", "Si2-alt.cell"
        )
        self.zns_band_cell = os.path.join(ilr_files("tests"), "data", "ZnS", "zns.cell")
        self.zns_singlepoint_cell = os.path.join(
            ilr_files("tests"), "data", "ZnS", "zns-sp.cell"
        )
        si_structure_file = os.path.join(ilr_files("tests"), "data", "Si", "Si8.json")
        self.si_structure = Structure.from_file(si_structure_file)

    def test_castep_cell_null_init(self):
        null_cell = CastepCell()
        self.assertEqual(null_cell.blocks, {})
        self.assertEqual(null_cell.tags, {})

        with self.assertRaises(ValueError):
            null_cell.structure

    def test_castep_cell_from_singlepoint_file(self):
        cc = CastepCell.from_file(self.zns_singlepoint_cell)
        self.assertEqual(
            set(cc.blocks.keys()),
            {"lattice_cart", "positions_abs", "species_pot"},
        )
        self.assertEqual(
            {k: v.value for k, v in cc.tags.items()},
            {
                "fix_all_cell": ["true"],
                "fix_all_ions": ["true"],
                "symmetry_generate": ["true"],
                "kpoint_mp_grid": ["4", "4", "4"],
                "snap_to_symmetry": ["true"],
            },
        )
        self.assertEqual(cc.blocks["species_pot"].values[1], ["S", "NCP"])
        self.assertEqual(cc.blocks["species_pot"].comments, ["", ""])

        structure = cc.structure
        self.assertIsInstance(structure, Structure)
        assert_array_almost_equal(
            structure.lattice.matrix,
            [[0.0, 2.71, 2.71], [2.71, 0.0, 2.71], [2.71, 2.71, 0.0]],
        )

    def test_castep_cell_from_structure(self):
        cell = CastepCell.from_structure(self.si_structure)
        self.assertEqual(
            cell.blocks["lattice_cart"].values,
            [
                ["ang"],
                ["5.43", "0.0", "0.0"],
                ["0.0", "5.43", "0.0"],
                ["0.0", "0.0", "5.43"],
            ],
        )
        self.assertEqual(
            cell.blocks["positions_frac"].values[0], ["Si", "0.0", "0.0", "0.0"]
        )
        self.assertEqual(
            cell.blocks["positions_frac"].values[7],
            ["Si", "0.75", "0.75", "0.25"],
        )


class CastepDosNiOTestCase(unittest.TestCase):
    def setUp(self):
        nio_files = {
            "bands_file": "NiO.bands",
            "pdos_file": "NiO.pdos_bin",
            "cell_file": "NiO.cell",
        }
        for key, value in nio_files.items():
            nio_files[key] = os.path.join(ilr_files("tests"), "data", "NiO", value)

        self.read_dos_kwargs = nio_files

        json_files = {"dos": "dos-pmg.json.gz", "pdos": "pdos-pmg.json.gz"}
        self.ref_data = {}

        for key, value in json_files.items():
            filename = os.path.join(ilr_files("tests"), "data", "NiO", value)

            with gzip.open(filename) as f:
                self.ref_data[key] = json.load(f, cls=MontyDecoder)

    def test_tdos_only(self):
        kwargs = self.read_dos_kwargs.copy()
        del kwargs["pdos_file"]
        del kwargs["cell_file"]

        tdos, pdos = read_dos(**kwargs)

        for spin, densities in self.ref_data["dos"].densities.items():
            assert_array_almost_equal(tdos.densities[spin], densities)
        self.assertFalse(pdos)

    def test_missing_cell_file(self):
        """Check error raised if PDOS given without .cell file"""
        kwargs = self.read_dos_kwargs.copy()
        del kwargs["cell_file"]

        with self.assertRaises(OSError):
            tdos, pdos = read_dos(**kwargs)

    def test_pdos(self):
        tdos, pdos = read_dos(**self.read_dos_kwargs)

        for spin, densities in self.ref_data["dos"].densities.items():
            assert_array_almost_equal(tdos.densities[spin], densities)

        for component, orbital_data in self.ref_data["pdos"].items():
            for orbital, dos in orbital_data.items():
                for spin, densities in dos.densities.items():
                    assert_array_almost_equal(
                        pdos[component][orbital].densities[spin], densities
                    )


class CastepBandStructureTestCaseNoSpin(unittest.TestCase):
    def setUp(self):
        self.si_bands = os.path.join(ilr_files("tests"), "data", "Si", "Si2.bands")
        self.si_cell = os.path.join(ilr_files("tests"), "data", "Si", "Si2.cell")
        self.si_cell_alt = os.path.join(
            ilr_files("tests"), "data", "Si", "Si2-alt.cell"
        )
        self.si_header_ref = os.path.join(
            ilr_files("tests"), "data", "Si", "Si2.bands_header.json"
        )

        self.ref_labels = {
            r"\Gamma": (0.0, 0.0, 0.0),
            "W": (0.5, 0.25, 0.75),
            "L": (0.5, 0.5, 0.5),
            "X": (0.5, 0.0, 0.5),
            "K": (0.375, 0.375, 0.75),
        }

    def test_castep_bands_read_header(self):
        header = read_bands_header(self.si_bands)
        with open(self.si_header_ref) as f:
            ref_header = json.load(f)
        self.assertEqual(header, ref_header)

    def test_castep_bands_read_eigenvalues(self):
        with open(self.si_header_ref) as f:
            ref_header = json.load(f)
        kpoints, weights, eigenvals = read_bands_eigenvalues(self.si_bands, ref_header)

        for i, k in enumerate([0.5, 0.36111111, 0.63888889]):
            self.assertAlmostEqual(kpoints[4][i], k)

        self.assertAlmostEqual(eigenvals[Spin.up][2, 4], 0.09500443 * _ry_to_ev * 2)

        for weight in weights:
            self.assertAlmostEqual(weight, 0.02272727)

    def test_castep_cell_read_labels(self):
        labels = labels_from_cell(self.si_cell)
        self.assertEqual(labels, self.ref_labels)

    def test_castep_cell_read_labels_alt_spelling(self):
        # CASTEP input is case-insensitive anf allows BS_KPOINTS_LIST
        # as well as BS_KPOINT_LIST                            ^
        labels = labels_from_cell(self.si_cell_alt)
        self.assertEqual(labels, self.ref_labels)


class CastepBandStructureTestCaseNickel(unittest.TestCase):
    def setUp(self):
        self.ni_cell = os.path.join(ilr_files("tests"), "data", "Ni", "ni-band.cell")

        self.ref_labels = {
            r"\Gamma": (0.0, 0.0, 0.0),
            "L": (0.5, 0.5, 0.5),
            "W": (0.5, 0.25, 0.75),
            "X": (0.5, 0.0, 0.5),
        }

    def test_castep_cell_read_labels_from_list(self):
        # The Si example uses handwritten .cell files in line-mode.
        # This example using a kgen-written .cell file with all k-points listed
        labels = labels_from_cell(self.ni_cell)
        self.assertEqual(labels, self.ref_labels)


class CastepBandStructureTestCaseWithSpin(unittest.TestCase):
    def setUp(self):
        self.fe_bands = os.path.join(ilr_files("tests"), "data", "Fe", "Fe.bands")

        self.fe_cell = os.path.join(ilr_files("tests"), "data", "Fe", "Fe.cell")
        self.fe_header_ref = os.path.join(
            ilr_files("tests"), "data", "Fe", "Fe.bands_header.json"
        )

    def test_castep_bands_read_header(self):
        header = read_bands_header(self.fe_bands)
        with open(self.fe_header_ref) as f:
            ref_header = json.load(f)
        self.assertEqual(header, ref_header)


class BandStructureTestCasePathBreak(unittest.TestCase):
    def setUp(self):
        self.pt_cell = os.path.join(ilr_files("tests"), "data", "Pt", "Pt.cell")
        self.ref_labels = {
            r"\Gamma": (0.0, 0.0, 0.0),
            "X": (0.5, 0.0, 0.5),
            "U": (0.625, 0.25, 0.625),
            "K": (0.375, 0.375, 0.75),
            "L": (0.5, 0.5, 0.5),
            "W": (0.5, 0.25, 0.75),
        }

    def test_castep_parse_break_in_k_path(self):
        labels = labels_from_cell(self.pt_cell)
        print(labels)
        self.assertEqual(labels, self.ref_labels)


# To generate reference file::
#
#   sumo-phonon-bandplot -f zns.phonon --units cm-1 --to-json zns_phonon.json
class CastepPhononTestCaseZincblende(unittest.TestCase):
    def setUp(self):
        self.zns_phonon = os.path.join(ilr_files("tests"), "data", "ZnS", "zns.phonon")
        self.zns_cell = os.path.join(ilr_files("tests"), "data", "ZnS", "zns.cell")
        self.zns_phonon_ref = os.path.join(
            ilr_files("tests"), "data", "ZnS", "zns_phonon.json"
        )

    def test_castep_phonon_read_bands(self):
        castep_phonon = CastepPhonon.from_file(self.zns_phonon)
        castep_phonon.set_labels_from_file(self.zns_cell)

        bs = castep_phonon.get_band_structure()
        bs_dict = bs.as_dict()

        with open(self.zns_phonon_ref) as f:
            ref_dict = json.load(f)

        for key in bs_dict.keys():
            self.assertIn(key, ref_dict.keys())

        # Comparing this data is quite brittle so compare the arrays when possible
        assert_array_almost_equal(
            bs_dict["lattice_rec"]["matrix"], ref_dict["lattice_rec"]["matrix"]
        )
        assert_array_almost_equal(bs_dict["qpoints"], ref_dict["qpoints"])
        self.assertEqual(bs_dict["labels_dict"], ref_dict["labels_dict"])
        self.assertEqual(bs_dict["structure"]["sites"], ref_dict["structure"]["sites"])
        assert_array_almost_equal(
            bs_dict["structure"]["lattice"]["matrix"],
            ref_dict["structure"]["lattice"]["matrix"],
        )

        assert_array_almost_equal(bs_dict["bands"], ref_dict["bands"])
        assert_array_almost_equal(
            bs_dict["eigendisplacements"]["real"],
            ref_dict["eigendisplacements"]["real"],
        )
        assert_array_almost_equal(
            bs_dict["eigendisplacements"]["imag"],
            ref_dict["eigendisplacements"]["imag"],
        )

        # May be an idea to check has_nac if that logic ever gets sorted out
