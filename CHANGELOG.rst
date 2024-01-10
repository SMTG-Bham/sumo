Change Log
==========

v2.3.8
------

Bugfixes:

- Fixed bug with sumo-phonon-bandplot and recent versions of pymatgen (@utf, #227)

v2.3.7
------

* Added symprec functionality to load_phonon in cli/phonon_bandplot.py by @badw in https://github.com/SMTG-Bham/sumo/pull/207
* Set band-edge marker size in style file by @ajjackson in https://github.com/SMTG-Bham/sumo/pull/193
* Fix failing tests and migrate to importlib by @oashour in https://github.com/SMTG-Bham/sumo/pull/216
* Update README.rst by @alexsquires in https://github.com/SMTG-Bham/sumo/pull/224
* Add support for band structure titles by @oashour in https://github.com/SMTG-Bham/sumo/pull/215
* Add ``**kwargs`` to ``draw_themed_line`` by @kavanase in https://github.com/SMTG-Bham/sumo/pull/212
* Kpoint degeneracy by @kbspooner in https://github.com/SMTG-Bham/sumo/pull/213

v2.3.6
------


Bugfixes:

- Allow "BREAK" keyword in CASTEP band paths - and ignore it. (@azanre, #201)
- Fix superscript in ``sumo-optplot`` (@utf,  #203)

v2.3.5
------

Fix projected band structures with two projections selected (@utf, #181).

v2.3.4
------

Support for new versions of `castepxbin`.

v2.3.3
------

Bugfix: ``sumo-bandplot --project --mode rgb`` now works with paths containing branches
with a single k-point.

v2.3.2
------

Bugfix: ``sumo-bandstats`` now reports the correct k-point indices for the VBM and CBM.


v2.3.1
------

You can now specify custom colours for RGB mode projected band structures. See the
documentation for more details. The colormath package has been added as a dependency.

v2.3.0
------

Major changes:

- DOS energies are no longer shifted by SIGMA in smeared calculations;
  this will cause the DOS line to extend slightly beyond the VBM, but
  will ensure peaks are in the right positions. (@yw-fang & @ajjackson)

New features:

- Fermi level may be chosen as energy zero in DOS or band structure (YWF, AJJ)
- Horizontal line may be requested at energy zero in DOS or band structure (YWF, AJJ)

Bugfixes

- VBM shift argument was not correctly passed in DOS plots (@pzarabadip)
- Scissor option with combined band and DOS plot did not scissor DOS (@utf)
- Spin option with combined band and DOS plot did not apply to the DOS (@utf)
- Fixed phonon band structures with non-analytic correction (@utf)

v2.2.5
------

Bugfix: Fixed colour parsing for BS + DOS plots. (@kavanse)

v2.2.4
------

Bugfix: ``sumo-kgen`` with ``--hybrid`` but without ``--folders`` now asks whether to
split up the k-points.

v2.2.3
------

Bugfix: Questaal band structure import was consistently raising unnecessary Exception

v2.2.2
------

Added Latimer-Munro high-symmetry k-point paths. (@kavanse)

v2.2.1
------

Fixed typo in ``sumo-bandplot``.

v2.2.0
------

- Support is added for orbital-projected DOS plots from CASTEP. (@zhubonan)

  - The binary parser is implemented and maintained in a separate
    library castepxbin. This is maintained by Bonan Zhu, available on
    PyPI and pinned to a specific version in the Sumo setup.py.

Enhancements:

- ``normalise`` option added to ``bandplot`` to control the normalisation of orbital
  projections. The default has been changed from ``select`` to ``all``, meaning that
  the size of projections is normalised against the sum of all other projections
  at that band and k-point. (@utf)

v2.1.1
------

Enhancements:

- Band structure grid lines can now be customised using a matplotlib
  style sheet.

v2.1.0
------

Sumo is now python 3.6+ only.

Additional bug fixes:

- Fix band structure interpolation with small branches (@kavanase)
- Update pymatgen version requirement.


v2.0.2
------

New testing and release framework.

v2.0.1
------

Bug fixes:

- Fixed support for pymatgen versions > 2020.10.9.1 (@utf)
- Fix yaml phonon-bandstructure plotting (@kavanase)


v2.0.0
------

New features:

- Support for CASTEP: (AJJ)

  - kgen: reciprocal-space path generation for electronic and phonon band-structures
  - bandplot: band structures (with or without spin-polarisation). Currently no
    support for element/orbital projected data (which would require a binary file parser).
  - phonon-bandplot: phonon band structures from .phonon files
  - dosplot: total-DOS plotting from eigenvalues. Again,
    projected-DOS plots are not currently available.

Bug fixes:

- Fix an oversight in the initial CASTEP/kgen implementation when the user provides a non-primitive cell as input.
- Python API fix for spin selection. (@kavanase)
- Fix phonon band structure line density selection. (@utf)

v1.4.0
------

New features:

- Plot single spin channel band structures. (@kavanase)
- Add scissor option to band plot. (@mkhorton)

Bug fixes:

- Fixed ytick labels for band + DOS plots. (@utf)
- Fix a bug when the y axis limit is outside the DOS range in band + DOS tapes. (@utf)

v1.3.0
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

v1.2.0
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

v1.1.3
------

Update Manifest.in

v1.1.2
------

Various bugfixes and enhancements:

- Fix manual k-point selection in kgen.
- Band indicies in bandstats now 1-based.
- Fix colour cycler issue in band structures with DOS.
- Allow overriding y-axis DOS ticks.
- Fermi level now set to 0 eV in dosplot .dat files (@frssp).
- Add ``--units`` option for phonon band structures (@ajjackson).
- Remove numbers from x-axis in band structures with DOS.

v1.1.1
------

Fix bug when installing from Pypi.


v1.1.0
------

Use matplotlib style sheets for styling plots (@ajjackson & @utf).
Enables plots to be customised based on user settings.

Various bug fixes:

- Fix bug when normalising DOS to Fermi level.
- Fix codacy style issues.
- Plotting style standardised across all plots.

v1.0.10
-------

Add option to align DOS to Fermi level (@shyamd)

Various bug fixes:

- Fix many typos.
- Updates to paper and documentation.

v1.0.9
------

``phonon-bandplot`` now supports combined DOS & band structure plots (Adam Jackson, Arthur Yaud).

Various bug fixes:

- Fix P centered trigonal k-point path.
- Fix ``--symprec`` behaviour in phonon-bandplot.
- Fix orbital projected band structures with branches (Adam Jackson).
- Fix reading Eg from spin-pol calculations (Adam Jackson).

v1.0.8
------

Enhancements by Adam Jackson:

- Add y-label and dos label options for DOS & band plots.
- Cache DOS colours for consistent plots.

Various bug fixes:

- Fixed gaussian broadening of DOS.
- Fixed ``--spg`` option in kgen and phonon-bandplot.
- Fixed default arguments for band structure + dos plotting.
- Added A centered orthorhombic lattice to ``BradCrackKpath``.

v1.0.7
------

Various bug fixes:

- Fixed density option in kgen.
- Fixed phonon-bandplot plotting limits.

v1.0.6
------

Move package data files.

v1.0.5
------

Minor bug fixes.

v1.0.4
------

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
