Change Log
==========

`[Unreleased] <https://github.com/smtg-ucl/sumo/compare/v1.3.0...HEAD>`_
------------------------------------------------------------------------

Major changes:

- DOS energies are no longer shifted by SIGMA in smeared calculations;
  this will cause the DOS line to extend slightly beyond the VBM, but
  will ensure peaks are in the right positions. (@yw-fang & @ajjackson)

New features:

- When plotting spin-polarised DOS, can choose to plot one channel only by specifying spin (@kavanase)
- Fermi level may be chosen as energy zero in DOS or band structure (YWF, AJJ)
- Horizontal line may be requested at energy zero in DOS or band structure (YWF, AJJ)


`[v1.3.0] <https://github.com/smtg-ucl/sumo/compare/v1.2.0...v1.3.0>`_ - 2019-10-20
------

This is the last supported version for Python 3.5, due to changes in pymatgen.

New features:

- Ability to plot multiple phonon band structures on top of each other. (AJJ)
- Added primitive-auto option to ``sumo-phonon-bandplot``. AJJ

Bug fixes:

- Added compatability with matplotlib 3.1. (AJJ)
- Use primitive cell when reading BORN. (AJJ)
- Set DOS cutoff when using ``--no-total``. (AJJ)
- Fix custom styling for phonon bandplotting from the CLI. (AJJ)
- Fix rare interpolation issues for projected band structure plots. (@utf)

`[v1.2.0] <https://github.com/smtg-ucl/sumo/compare/v1.1.3...v1.2.0>`_ - 2019-06-07
------

This is the most contributers to a release so far!

- Bug fixes
  - Fix error in P monoclinic (*b*-unique) "Bradcrack" high-symmetry path (C. N. Savory)
  - Fix appearance of y-axis formatter for optics plots (E. Rubinstein & Adam J. Jackson)
  - Prevent an error when requesting DOS subplots with no total DOS (Z. Xing)
  - Fix missing f0 orbitals in orbital projected DOS plots (@utf)
  - Update phonon-bandplot to use latest phonopy API (@utf)

- New features

  - Additional properties from dielectric function (AJJ & K. T. Butler)

    - any combination of absorption, loss, dielectric and complex refractive index components can be requested as a set of subplots

  - Allow full 3x3 supercell matrix to be specified for phonon band structures (AJJ)

  - Band structure label manipulation with '@' (AJJ)

    - place @ before a label to make it invisible in plot
    - place @ at end to make unique point that avoids confusing pymatgen; the label will be displayed without any trailing @ characters

  - Aspect ratio control for band structures (A. M. Ganose)

- New interfaces

  - Questaal is now supported. Pretty much everything works except
    orbital-decomposed band structures and phonons. (AJJ)

    - Generate a *syml.ext* band path file using **sumo-kgen** reading from a
      site.ext or init.ext file. (The site file is generally the correct
      choice.)

    - Plot electronic band structure generated with **lmf** using
      **sumo-bandplot** reading from *bnds.ext* and *syml.ext* files.

    - Plot a total DOS from **lmf** with **sumo-dosplot** reading *ext.dos*

    - Plot a PDOS from **lmf** by moving the total dos to *tdos.ext*
      and using Questaal tools to generate a *dos.ext* with orbital
      information before running **sumo-dosplot**.

    - Plot optical properties with **sumo-optplot** from dielectric
      function written by **lmf** (*opt.ext*) or **bethesalpeter**
      (*ext.eps_BSE*). Optical spectra from multiple sources
      (e.g. VASP and bethesalpeter) may be plotted alongside one
      another.

`[v1.1.3] <https://github.com/smtg-ucl/sumo/compare/v1.1.2...v1.1.3>`_ - 2018-12-24
-----------------------------------------------------------------------------------

Update Manifest.in

`[v1.1.2] <https://github.com/smtg-ucl/sumo/compare/v1.1.1...v1.1.2>`_ - 2018-11-30
-----------------------------------------------------------------------------------

Various bugfixes and enhancements:

- Fix manual k-point selection in kgen.
- Band indicies in bandstats now 1-based.
- Fix colour cycler issue in band structures with DOS.
- Allow overriding y-axis DOS ticks.
- Fermi level now set to 0 eV in dosplot .dat files (@frssp).
- Add ``--units`` option for phonon band structures (@ajjackson).
- Remove numbers from x-axis in band structures with DOS.

`[v1.1.1] <https://github.com/smtg-ucl/sumo/compare/v1.1.0...v1.1.1>`_ - 2018-08-15
-----------------------------------------------------------------------------------

Fix bug when installing from Pypi.


`[v1.1.0] <https://github.com/smtg-ucl/sumo/compare/v1.0.10...v1.1.0>`_ - 2018-08-10
------------------------------------------------------------------------------------

Use matplotlib style sheets for styling plots (@ajjackson & @utf).
Enables plots to be customised based on user settings.

Various bug fixes:

- Fix bug when normalising DOS to Fermi level.
- Fix codacy style issues.
- Plotting style standardised across all plots.

`[v1.0.10] <https://github.com/smtg-ucl/sumo/compare/v1.0.9...v1.0.10>`_ - 2018-08-06
-------------------------------------------------------------------------------------

Add option to align DOS to Fermi level (@shyamd)

Various bug fixes:

- Fix many typos.
- Updates to paper and documentation.

`[v1.0.9] <https://github.com/smtg-ucl/sumo/compare/v1.0.8...v1.0.9>`_ - 2018-06-19
-----------------------------------------------------------------------------------

``phonon-bandplot`` now supports combined DOS & band structure plots (Adam Jackson, Arthur Yaud).

Various bug fixes:

- Fix P centered trigonal k-point path.
- Fix ``--symprec`` behaviour in phonon-bandplot.
- Fix orbital projected band structures with branches (Adam Jackson).
- Fix reading Eg from spin-pol calculations (Adam Jackson).

`[v1.0.8] <https://github.com/smtg-ucl/sumo/compare/v1.0.7...v1.0.8>`_ - 2018-05-9
----------------------------------------------------------------------------------

Enhancements by Adam Jackson:

- Add y-label and dos label options for DOS & band plots.
- Cache DOS colours for consistent plots.

Various bug fixes:

- Fixed gaussian broadening of DOS.
- Fixed ``--spg`` option in kgen and phonon-bandplot.
- Fixed default arguments for band structure + dos plotting.
- Added A centered orthorhombic lattice to ``BradCrackKpath``.

`[v1.0.7] <https://github.com/smtg-ucl/sumo/compare/v1.0.6...v1.0.7>`_ - 2018-04-24
-----------------------------------------------------------------------------------

Various bug fixes:

- Fixed density option in kgen.
- Fixed phonon-bandplot plotting limits.

`[v1.0.6] <https://github.com/smtg-ucl/sumo/compare/v1.0.5...v1.0.6>`_ - 2018-04-18
-----------------------------------------------------------------------------------

Move package data files.

`[v1.0.5] <https://github.com/smtg-ucl/sumo/compare/v1.0.4...v1.0.5>`_ - 2018-04-17
-----------------------------------------------------------------------------------

Minor bug fixes.

`[v1.0.4] <https://github.com/smtg-ucl/sumo/compare/v1.0.0...v1.0.4>`_ - 2018-04-16
-----------------------------------------------------------------------------------

Minor changes to Pypi config.

v1.0.0
------

Added
~~~~~

- Script files:

  - sumo-kgen
  - sumo-dosplot
  - sumo-bandplot
  - sumo-bandstats
  - sumo-optplot
  - sumo-phonon-bandplot

