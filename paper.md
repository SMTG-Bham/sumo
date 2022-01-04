---
title: 'sumo: Command-line tools for plotting and analysis of periodic *ab initio* calculations'
tags:
  - plotting
  - ab initio
  - density of states
  - band structure
  - optical absorption
  - effective mass
  - phonons
  - vasp
authors:
 - name: Alex M Ganose
   orcid: 0000-0002-4486-3321
   affiliation: "1, 2, 3"
 - name: Adam J Jackson
   orcid: 0000-0001-5272-6530
   affiliation: 1
 - name: David O Scanlon
   orcid: 0000-0001-9174-8601
   affiliation: "1, 2, 3"
   email: d.scanlon@ucl.ac.uk
affiliations:
 - name: Dept of Chemistry, University College London, 20 Gordon Street, London WC1H 0AJ, UK
   index: 1
 - name: Diamond Light Source Ltd., Diamond House, Harwell Science and Innovation Campus, Didcot, Oxfordshire OX11 0DE, UK
   index: 2
 - name: Thomas Young Centre, University College London, Gower Street, London WC1E 6BT, UK
   index: 3
date: April 2018
bibliography: paper.bib
---

*Ab initio* electronic structure modelling is capable of providing an
insight into the fundamental properties of solid-state materials, at a resolution
beyond that of experimental techniques. The optoelectronic properties
of a compound are analysed through several key descriptions, including:
density-of-states distributions, which provide information on the
orbital character of bonding; band structure diagrams, which indicate
carrier transport properties; and optical absorption spectra, which are
used to assess the wavelengths of light a material will transmit or
absorb. An understanding of these fundamental properties is crucial
when selecting or optimising materials for particular applications,
including photovoltaics [@solar], transparent conductors [@TCO], and
thermoelectrics [@thermoelectrics].

Most common *ab initio* calculation software for analysing crystalline
materials with periodic boundary condictions, such as Vienna *ab
initio* Simulation Package (VASP) [@vasp] and Quantum Espresso
[@QEcode], write raw data which require post-processing to plot or
convert into a human-readable format. Several packages exist that
facilitate the creation and plotting of such diagrams. Python
libraries, such as Python Materials Genomics (pymatgen) [@pymatgen]
and Atomic Simulation Environment (ase) [@ase], provide powerful
interfaces for plotting and data analysis but require the user to be
proficient in Python to use effectively. Conversely, programs which
provide a graphical user interface, such as p4vasp [@p4vasp] and
XCrySDen [@xcrysden], are easy to use but are not conducive to working
on the command line. The purpose of this package is to provide an
intermediate solution that is trivial to use but still provides the
flexibility needed for a broad range of analysis modes.


# `sumo`

`sumo` is a set of command-line tools for publication-ready plotting
and analysis of *ab initio* calculation data for solid-state materials.
The code includes a
fully-documented Python module, upon which the command-line
scripts are built. `sumo` currently only supports VASP, however,
extending the code to other solid-state *ab initio* calculators is planned for future
releases. The code relies on several open-source Python packages for
common tasks, including pymatgen for data loading [@pymatgen], spglib
for symmetry analysis [@spglib], and matplotlib for plotting
[@matplotlib].

The main plotting functionality of `sumo` includes density of states
plots, electronic and phonon band structure diagrams, and optical
absorption spectra (Figure 1). The code has been designed to allow for
significant customisation of plots, including the ability to produce
projected density of states and orbital resolved band structures. The
code additionally supplies a tool for generating k-point paths along
high-symmetry directions in the Brillouin zone, with the ability to
write the necessary input files required to perform the
calculations in VASP. Crucially, this tool allows a single band
structure plot to be split into several *ab initio* calculations,
as is essential when dealing with large materials or restrictive batch
systems. Lastly, a script is provided to extract information from
semiconductor band structures, including direct and indirect band gaps,
band edge locations, and parabolic and non-parabolic effective masses.

![Diagrams produced by `sumo`. a) Density of states, b) projected band
structure, and c) optical absorption spectra.](docs/src/figures/sumo_plots.pdf)


# Acknowledgements

DOS acknowledges support from the EPSRC (EP/N01572X/1). DOS
acknowledges support from the European Research Council, ERC (grant no.
758345). DOS acknowledges membership of the Materials Design Network.
AMG acknowledges Diamond Light Source for the co-sponsorship of a
studentship on the EPSRC Centre for Doctoral Training in Molecular
Modelling and Materials Science (EP/L015862/1).

We acknowledge useful discussions with Zhenyu Wang, Benjamin Morgan,
and Jonathan Skelton. Feature requests and user testing came from
Benjamin Williamson, Christopher Savory and James Pegg.


# References
