ScriptsForVASP
==============
Making life easier using scripting languages (Bash and Python) to facilitate multiple VASP simulation job preparation, submission and analysis.

Introduction
--------------

There are multiple ways to carry out VASP calculations.

The somewhat "primitive" way is creating a directory, writing up or copying from somewhere else the input files, mostly INCAR, KPOINTS, POSCAR, POTCAR, and a job submission file depending on the supercomputer one is using, then submiting the job and waiting for the output files showing up in the same directory. 

But most of the time, one single run is far from enough to extract the information one wishes to obtain. Rather, a series of runs with one or more varying parameters are needed, a post-analysis should be done, and a cumulative but terse output in the form of a file or direct screen output should be provided to the user. Examples are __cut-off energy and k-points convergence tests__, __lattice constant (equilibrium volume) and bulk modulus test__, __elastic constant test__, and some __exotic tests__ that require a __gradual change__ somewhere in the __POSCAR__, such as a __rotation of several chosen ions__. No one should expect a user to have the patience to change the varying parameter of every single run in the series and submit them one by one, extract information from each run by hand, finally fit the curve, conclude the coefficients and preferably plot them in a saved figure. In order to tackle this problem, there are two automated ways.

The __first__ is wrapping up the job preparation, execution and post-analysis trilogy in an installable package, written in and for a high level scripting language, such as Python. A famous example is [ASE] (https://wiki.fysik.dtu.dk/ase/) (Atomic Simulation Environment). One needs to install the Python package and read the documentation on the location and functionality of their defined classes and methods, write up a Python script for a series of runs, oftentimes along with the post-analysis part and feed it to the batch job submission script. It has a lot of advantages, especially if one deals with different kinds of _ab initio_ calculators and very routine work that all of them share, because it hides the basic file and directory manipulation, which are different for different calculators, and let users work with its defined Python classes and methods.

The __second__ is providing a collection of scripts, written in a shell scripting language, like Bash, and accompanying them with scripts written in a programming scripting language, like Python, in the same directory, with careful use of file names and input arguments, and letting users add this single directory into their PATH environmental variable to have direct access. The scripts will be used as executables. This project aims to strengthen such an approach. This approach is less systematic than the former, but it also has some obvious advantages comparitively.

* It provides more flexibility in terms of file and directory manipulation, because Bash, as a shell language, deals with them more natively than Python, if one doesn't bother to use more Python packages in order to replace Bash. The directory that ASE runs in will not hold the full output of a series of runs, but only the last one, in an __equilibrium volume and bulk modulus__ routine. But with Bash one can conveniently create subdirectories to hold the input and output of each run in an automated way, thus keeping records of everything for future use.

* With all of the output files present, one can do far more checks and data collection with some knowledge of command line tools such as __cat__, __grep__, __sed__, __awk__, and __find__. These tools are considered more universal than the layered classes and methods defined in ASE Python package. The executable scripts in a single directory are grouped more straighforwardly, showing the users only task related commands with simple keywords. As far as my own experience, VASP is already a delicate machine that requires enough memorization of its tags, and through reading the docs on them written by the VASP creators, I automatically get the idea of where to look for the desired information in the output files. There is really less motivation to learn a completely new Pythonic scheme. Plus, if the physics I'm looking at needs custom expression and control than standard methods provided by a Python package, I'm more likely to express it myself with __grep__, __sed__ and __awk__, than wait for the developers of ASE to come up with a new method and a new code block somewhere in their documentation. At least, at this stage in my opinion, there are more physics than the convenient expression of Python classes and methods.

* One can freely separate processes of the input files preparation, job submission and post-analysis, just as implemented in this project by __Prepare.sh__, __Fire.sh__, and __Display.sh__, or combine the first two for daily use when feeling secure as __Prep-fire.sh__. For a single run, one can just use __Fast-prep.sh__. The overall directory structure will be very neat. This opens a door widely to runs that depend on the output of previous runs. It's just several lines of commands away with all the convenient tools in the linux shell.

* If one writes up a Python script describing a series of runs with ASE, and submit the whole thing, and it reaches the walltime, or one or more of the runs crash, or the script is faulty somewhere, which are scenarios that always happen, the script cannot be finished, and it would be difficult to pick up what has been terminated unexpectedly. The scheme is somewhat fragile. With the existence of multiple subdirectories, however, one can give each run in a series an independent job submission script, and submit them selectively. If one run reaches the walltime or simply fails, one can easily identify it, change the input parameter and submit it again, and other runs won't be affected. If there are mistakes in the job preparation scripts, one will discover them before submission, because chances are that those mistakes are directly output on the screen before one fire up the submission scripts.


* To save more space, one can choose not to output, or just delete the WAVECAR, CHG, and CHGCAR, which are disk space eaters.

Based on above reasons, which of course will be less persuasive for an experience ASE user who is already familiar with its Pythonic scheme and capable of freely manipulating directories anyway, or if ASE improves its management of directories, I present the following scripts that are used by me on a daily basis for dealing with routines described in the beginning of this section, and scripts offering more functionalities, including the usual need for __density of states plots__, __band structure plots__, __Bader charge transfer analysis__, etc. The list will go on and can be mixed with other people's scripts conveniently.

Due to the low level of "abstraction", it is assumed that the user of this script collection is familiar with simple Bash syntax, command line tools such as __cat__, __grep__, __sed__, __awk__, __find__ and some Python. Flexibility comes with some price, so one may need to look at the codes and change some lines. But since one doesn't have to be forced to face a formidable Python package with intricate definitions, it will be much easier. I'll try to turn the hard coded parts to automated detection and sweet input command line arguments, and make appropriate comments on the routines as I progress. This is still an ongoing project, but you are welcome to fork it and modify the scripts to your own need.

Scripts structure
---------------------

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
