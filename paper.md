---
title: 'sumo: Command-line plotting tools for ab initio calculations'
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
   affiliation: 1,2,3
 - name: Adam J Jackson
   orcid: 0000-0001-5272-6530
   affiliation: 1
 - name: David O Scanlon
   orcid: 0000-0001-9174-8601
   affiliation: 1,2,3
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

Ab initio electronic structure modelling is capable of providing an insight into the
fundamental properties of materials, at a resultion beyond that of experimental techinques.
The optoelectronic properties of a compound can be described through several key
descriptions, including: density of states diagrams, which provide information on the
orbital character of frontier bonding; band structures, which afford an indication of
carrier transport properties; and optical absorption spectra, which can be used to
assess the wavelengths of light a material will transmit or absorb.
An understanding of these fundamental properties is crucial when selecting or optimising
materials for particular applications, including photovoltaics, transparent conductors, and
thermoelectrics.

Most common ab initio calculation software, such as Vienna ab initio Simulation Package (VASP)
and Quantum Espresso, write raw data which require post-processing steps to plot or convert
into a human-readable format.
Several packages exist that facilitate the creation and plotting of such diagrams.
Python libraries, such as Python Materials Genomics (pymatgen) and Atomic Simulation
Environment (ase), provide powerful interfaces for plotting and data analysis but
require the user to be profficient in Python to use effectively.
Conversely, programs which provide a graphical user interface, such as p4vasp and XCrySDen,
are easy to use but are not conducive to working on the command line.

# `sumo`

`sumo` is a set of command line scripts and a corresponding Python API, for publication-ready
plotting and analysis of ab initio calculation data.

- publication ready diagrams
- customisable
- quick and easy to use
- python api available
- list of commands available
- example of band structure, dos and optics plots
