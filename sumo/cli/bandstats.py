# Copyright (c) Scanlon Materials Theory Group
# Distributed under the terms of the MIT License.

"""
A script to calculate effective masses and band degeneracies.

TODO:
 - Would be good to reimplement get_vbm and get_cbm methods, with the ability
   to give degeneracies + other band edges within a certain tolerance.
 - Sample custom k-point, band, spin combinations.
"""

import argparse
import logging
import sys

from pymatgen.electronic_structure.core import Spin
from pymatgen.io.vasp.outputs import BSVasprun

from sumo.cli.bandplot import find_vasprun_files
from sumo.electronic_structure.bandstructure import get_reconstructed_band_structure
from sumo.electronic_structure.effective_mass import (
    fit_effective_mass,
    get_fitting_data,
)

__author__ = "Alex Ganose"
__version__ = "1.0"
__maintainer__ = "Alex Ganose"
__email__ = "alexganose@googlemail.com"
__date__ = "March 4, 2018"


kpt_str = "[{k[0]:.2f}, {k[1]:.2f}, {k[2]:.2f}]"


def bandstats(
    filenames=None,
    num_sample_points=3,
    temperature=None,
    degeneracy_tol=1e-4,
    parabolic=True,
):
    """Calculate the effective masses of the bands of a semiconductor.

    Args:
        filenames (:obj:`str` or :obj:`list`, optional): Path to vasprun.xml
            or vasprun.xml.gz file. If no filenames are provided, the code
            will search for vasprun.xml or vasprun.xml.gz files in folders
            named 'split-0*'. Failing that, the code will look for a vasprun in
            the current directory. If a :obj:`list` of vasprun files is
            provided, these will be combined into a single band structure.
        num_sample_points (:obj:`int`, optional): Number of k-points to sample
            when fitting the effective masses.
        temperature (:obj:`int`, optional): Find band edges within kB * T of
            the valence band maximum and conduction band minimum. Not currently
            implemented.
        degeneracy_tol (:obj:`float`, optional): Tolerance for determining the
            degeneracy of the valence band maximum and conduction band minimum.
        parabolic (:obj:`bool`, optional): Use a parabolic fit of the band
            edges. If ``False`` then nonparabolic fitting will be attempted.
            Defaults to ``True``.

    Returns:
        dict: The hole and electron effective masses. Formatted as a
        :obj:`dict` with keys: ``'hole_data'`` and ``'electron_data'``. The
        data is a :obj:`list` of :obj:`dict` with the keys:

        'effective_mass' (:obj:`float`)
            The effective mass in units of electron rest mass, :math:`m_0`.

        'energies' (:obj:`numpy.ndarray`)
            Band eigenvalues in eV.

        'distances' (:obj:`numpy.ndarray`)
            Distances of the k-points in reciprocal space.

        'band_id' (:obj:`int`)
            The index of the band,

        'spin' (:obj:`~pymatgen.electronic_structure.core.Spin`)
            The spin channel

        'start_kpoint' (:obj:`int`)
            The index of the k-point at which the band extrema occurs

        'end_kpoint' (:obj:`int`)
    """
    if not filenames:
        filenames = find_vasprun_files()
    elif isinstance(filenames, str):
        filenames = [filenames]

    bandstructures = []
    for vr_file in filenames:
        vr = BSVasprun(vr_file, parse_projected_eigen=False)
        bs = vr.get_band_structure(line_mode=True, efermi="smart")
        bandstructures.append(bs)
    bs, kpt_mapping = get_reconstructed_band_structure(
        bandstructures, force_kpath_branches=True, return_forced_branch_kpt_map=True
    )

    if bs.is_metal():
        logging.error("ERROR: System is metallic!")
        sys.exit()

    _log_band_gap_information(bs, kpt_mapping=kpt_mapping)

    vbm_data = bs.get_vbm()
    cbm_data = bs.get_cbm()

    logging.info("\nValence band maximum:")
    _log_band_edge_information(bs, vbm_data, kpt_mapping=kpt_mapping)

    logging.info("\nConduction band minimum:")
    _log_band_edge_information(bs, cbm_data, kpt_mapping=kpt_mapping)

    if parabolic:
        logging.info("\nUsing parabolic fitting of the band edges")
    else:
        logging.info("\nUsing nonparabolic fitting of the band edges")

    if temperature:
        logging.error("ERROR: This feature is not yet supported!")

    else:
        # Work out where the hole and electron band edges are.
        # Fortunately, pymatgen does this for us. Points at which to calculate
        # the effective mass are identified as a tuple of:
        # (spin, band_index, kpoint_index)
        hole_extrema = []
        for spin, bands in vbm_data["band_index"].items():
            hole_extrema.extend(
                [
                    (spin, band, kpoint)
                    for band in bands
                    for kpoint in vbm_data["kpoint_index"]
                ]
            )

        elec_extrema = []
        for spin, bands in cbm_data["band_index"].items():
            elec_extrema.extend(
                [
                    (spin, band, kpoint)
                    for band in bands
                    for kpoint in cbm_data["kpoint_index"]
                ]
            )

        # extract the data we need for fitting from the band structure
        hole_data = []
        for extrema in hole_extrema:
            hole_data.extend(
                get_fitting_data(bs, *extrema, num_sample_points=num_sample_points)
            )

        elec_data = []
        for extrema in elec_extrema:
            elec_data.extend(
                get_fitting_data(bs, *extrema, num_sample_points=num_sample_points)
            )

    # calculate the effective masses and log the information
    logging.info("\nHole effective masses:")
    for data in hole_data:
        eff_mass = fit_effective_mass(
            data["distances"], data["energies"], parabolic=parabolic
        )
        data["effective_mass"] = eff_mass
        _log_effective_mass_data(data, bs.is_spin_polarized, mass_type="m_h")

    logging.info("\nElectron effective masses:")
    for data in elec_data:
        eff_mass = fit_effective_mass(
            data["distances"], data["energies"], parabolic=parabolic
        )
        data["effective_mass"] = eff_mass
        _log_effective_mass_data(data, bs.is_spin_polarized)

    return {"hole_data": hole_data, "electron_data": elec_data}


def _log_band_gap_information(bs, kpt_mapping=None):
    """Log data about the direct and indirect band gaps.

    Args:
        bs (:obj:`~pymatgen.electronic_structure.bandstructure.BandStructureSymmLine`):
        kpt_mapping (:obj:`dict`, optional): A mapping of k-point indicies from the
            band structure with forced branches to the original band structure.

    """
    bg_data = bs.get_band_gap()
    if not bg_data["direct"]:
        logging.info(f"Indirect band gap: {bg_data['energy']:.3f} eV")

    direct_data = bs.get_direct_band_gap_dict()
    if bs.is_spin_polarized:
        direct_bg = min(spin_data["value"] for spin_data in direct_data.values())
        logging.info(f"Direct band gap: {direct_bg:.3f} eV")

        for spin, spin_data in direct_data.items():
            direct_kindex = spin_data["kpoint_index"]
            direct_kpoint = bs.kpoints[direct_kindex].frac_coords
            direct_kpoint = kpt_str.format(k=direct_kpoint)
            eq_kpoints = bs.get_equivalent_kpoints(direct_kindex)
            eq_kpoints = _map_kpoints(eq_kpoints, kpt_mapping)
            k_indices = ", ".join(map(str, eq_kpoints))

            # add 1 to band indices to be consistent with VASP band numbers.
            b_indices = ", ".join([str(i + 1) for i in spin_data["band_indices"]])

            logging.info(f"  {spin.name.capitalize()}:")
            logging.info(f"    k-point: {direct_kpoint}")
            logging.info(f"    k-point indices: {k_indices}")
            logging.info(f"    Band indices: {b_indices}")

    else:
        direct_bg = direct_data[Spin.up]["value"]
        logging.info(f"Direct band gap: {direct_bg:.3f} eV")

        direct_kindex = direct_data[Spin.up]["kpoint_index"]
        direct_kpoint = kpt_str.format(k=bs.kpoints[direct_kindex].frac_coords)
        eq_kpoints = bs.get_equivalent_kpoints(direct_kindex)
        eq_kpoints = _map_kpoints(eq_kpoints, kpt_mapping)
        k_indices = ", ".join(map(str, eq_kpoints))
        b_indices = ", ".join(
            [str(i + 1) for i in direct_data[Spin.up]["band_indices"]]
        )

        logging.info(f"  k-point: {direct_kpoint}")
        logging.info(f"  k-point indices: {k_indices}")
        logging.info(f"  Band indices: {b_indices}")


def _log_band_edge_information(bs, edge_data, kpt_mapping=None):
    """Log data about the valence band maximum or conduction band minimum.

    Args:
        bs (:obj:`~pymatgen.electronic_structure.bandstructure.BandStructureSymmLine`):
            The band structure.
        edge_data (dict): The :obj:`dict` from ``bs.get_vbm()`` or
            ``bs.get_cbm()``
        kpt_mapping (:obj:`dict`, optional): A mapping of k-point indicies from the
            band structure with forced branches to the original band structure.
    """
    if bs.is_spin_polarized:
        spins = edge_data["band_index"].keys()
        b_indices = [
            ", ".join([str(i + 1) for i in edge_data["band_index"][spin]])
            + f"({spin.name.capitalize()})"
            for spin in spins
        ]
        b_indices = ", ".join(b_indices)
    else:
        b_indices = ", ".join([str(i + 1) for i in edge_data["band_index"][Spin.up]])

    kpoint = edge_data["kpoint"]
    kpoint_str = kpt_str.format(k=kpoint.frac_coords)
    k_indices = ", ".join(
        map(str, _map_kpoints(edge_data["kpoint_index"], kpt_mapping))
    )
    k_degen = bs.get_kpoint_degeneracy(kpoint=kpoint.frac_coords)

    if kpoint.label:
        k_loc = kpoint.label
    else:
        branch = bs.get_branch(edge_data["kpoint_index"][0])[0]
        k_loc = f"between {branch['name']}"

    logging.info(f"  Energy: {edge_data['energy']:.3f} eV")
    logging.info(f"  k-point: {kpoint_str}")
    logging.info(f"  k-point location: {k_loc}")
    logging.info(f"  k-point indices: {k_indices}")
    logging.info(f"  k-point degeneracy: {k_degen}")
    logging.info(f"  Band indices: {b_indices}")


def _log_effective_mass_data(data, is_spin_polarized, mass_type="m_e"):
    """Log data about the effective masses and their directions.

    Args:
        data (dict): The effective mass data. Formatted as a :obj:`dict` with
            the keys:

            'effective_mass' (:obj:`float`)
                The effective mass in units of electron rest mass, :math:`m_0`.

            'energies' (:obj:`numpy.ndarray`)
                Band eigenvalues in eV.

            'band_id' (:obj:`int`)
                The index of the band,

            'spin' (:obj:`~pymatgen.electronic_structure.core.Spin`)
                The spin channel

            'start_kpoint' (:obj:`int`)
                The index of the k-point at which the band extrema occurs

            'end_kpoint' (:obj:`int`)
                The k-point towards which the data has been sampled.

        is_spin_polarized (bool): Whether the system is spin polarized.
    """
    s = f" ({data['spin'].name})" if is_spin_polarized else ""

    # add 1 to band id to be consistent with VASP
    band_str = f"band {data['band_id'] + 1}{s}"

    start_kpoint = data["start_kpoint"]
    end_kpoint = data["end_kpoint"]
    eff_mass = data["effective_mass"]

    kpoint_str = kpt_str.format(k=start_kpoint.frac_coords)
    if start_kpoint.label:
        kpoint_str += f" ({start_kpoint.label})"
    kpoint_str += " -> "
    kpoint_str += kpt_str.format(k=end_kpoint.frac_coords)
    if end_kpoint.label:
        kpoint_str += f" ({end_kpoint.label})"

    logging.info(f"  {mass_type}: {eff_mass:.3f} | {band_str} | {kpoint_str}")


def _map_kpoints(kpt_idxs, kpt_mapping):
    """Map k-point indices to the original band structure."""
    if not kpt_mapping:
        return kpt_idxs
    return sorted(set([kpt_mapping.get(k, k) for k in kpt_idxs]))


def _get_parser():
    parser = argparse.ArgumentParser(
        description="""
    bandstats provides information on the band gap and effective
    masses of semiconductors.""",
        epilog=f"""
    Author: {__author__}
    Version: {__version__}
    Last updated: {__date__}""",
    )

    parser.add_argument(
        "-f",
        "--filenames",
        default=None,
        nargs="+",
        metavar="F",
        help="one or more vasprun.xml files to plot",
    )
    parser.add_argument(
        "-n",
        "--nonparabolic",
        default=True,
        action="store_false",
        help=("use a nonparabolic model to fit the effective masses"),
    )
    parser.add_argument(
        "-s",
        "--sample-points",
        default=3,
        type=int,
        dest="sample_points",
        metavar="N",
        help="number of k-points to sample in fitting",
    )
    return parser


def main():
    args = _get_parser().parse_args()
    logging.basicConfig(
        filename="sumo-bandstats.log",
        level=logging.INFO,
        filemode="w",
        format="%(message)s",
    )
    console = logging.StreamHandler()
    logging.info(" ".join(sys.argv[:]))
    logging.getLogger("").addHandler(console)

    bandstats(
        filenames=args.filenames,
        num_sample_points=args.sample_points,
        parabolic=args.nonparabolic,
    )


if __name__ == "__main__":
    main()
