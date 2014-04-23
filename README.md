ScriptsForVASP
==============
Making life easier using scripting languages (Bash and Python) to facilitate multiple VASP simulation job preparation, submission and analysis.

This project has a nicely formatted [document] (http://scriptsforvasp.readthedocs.org/). About anything below if you wish to read more, please have a look there.

Introduction
------------

There are multiple ways to carry out VASP calculations. Contrary to the Pythonic scheme the famous [ASE] (https://wiki.fysik.dtu.dk/ase/) uses, This project does things in an old-fashioned way, providing a collection of scripts, written in Bash, and accompanying them with scripts written in Python, in the same directory, with careful use of file names and input arguments, and letting users add this single directory into their PATH environmental variable to have direct access. The scripts will be used as executables. This approach is less systematic than the one ASE uses, but it is more file and directory aware, due to their somewhat natural affinity to Bash - a shell scripting language. A detailed analysis can be found [here] (http://scriptsforvasp.readthedocs.org/en/latest/introduction.html).

Scripts structure
-----------------

### Scripts that are of basic level

* `Fast-prep.sh XXX` - a fast way to create a directory XXX with necessary input files
* `Prepare.sh TEST_TYPE [OPTIONS]` - to create a series of runs, performing encut test / kpoints test / lattice constant test / elastic constant test / ...
* `Fire.sh TEST_TYPE [OPTIONS]` - to submit a series of jobs described above
* `Display.sh TEST_TYPE [OPTIONS]` - to do the post-analysis of a series of runs described above and determine physical quantities
* `Prep-fire.sh TEST_TYPE [OPTIONS]` - a wrapper on Prepare.sh and Fire.sh

### Scripts that are built upon the basic level

* `Elastic.sh prep-fire / disp-solve / solve` - to do a full set of independent elastic constant determination runs, built on `Prepare.sh`, `Fire.sh` and `Display.sh`
* `Dos-bs.sh scrun / dosrun / bsrun / plot-tdos / plot-ldos / plot-bs` - to carry out self-consistent runs, density of states runs, band structure runs, to plot the total and local density of states, and to plot the band structure
* `Bader.sh prerun / bader` - to perform Bader charge transfer analysis

### Scripts that are for multiple species, each requiring runs above

* `ElementRun.sh` - an interactive script to execute the above routines for multiple species

### Scripts that are in general useful and convenient

* `Cellinfo.sh [rwigs]` - to show the important info from OUTCAR (rwigs is optional, and displays the RWIGS values for different species to 100 % fill the cell)
* `Delfiles.sh [XXX]` - to delete the output files for failed runs (XXX is optional if one stays in the job directory)
* `Re-relax.sh [XXX]` - to copy CONTCAR to POSCAR and start over relaxation (XXX is optional if one stays in the job directory)
* `Swap-files-name.sh AAA BBB` - to exchange the names of two files AAA and BBB

### Python scripts that support the Bash scripst above

* `_Prepare-strain.py` - to add strain to the cell lattice vectors in POSCAR for a series of runs, determining the value of a combination of elastic constants
* `_Ions-rotator.py` - to rotate ions in selected lines in POSCAR at a certain point, around a certain axis, to a certain angle
* `_Display-fit.py` - to perform different fitting techniques, including Birch-Murnaghan equation of state and polynomial
* `_Elastic-solver.py` - to solve the independent elastic constants from the same number of combination runs with linear algebra
* `_Plot-tdos.py` - to plot total density of states, spin states are supported (above or below x axis)
* `_Plot-ldos.py` - to plot local projected density of states of two ions (above or below x axis), spin states are supported (in different files)
* `_Plot-bs.py` - to plot the band structure along certain chosen points in the 1st BZ
* `_Cellinfo-solver.py` - to obtain the RWIGS values for different species to 100 % fill the cell

__Non-job-related__ contains some interesting physics or computational experiments, and __Archived__ contains scripts that are obsolete or less used, but may contain jems.

General routines
---------------------
    HAVEN'T FINALIZED YET BELOW
