ScriptsForVASP
==============
Making life easier using scripting languages (Bash and Python) to facilitate multiple VASP simulation job preparation, submission and analysis.

Introduction
--------------

There are multiple ways to carry out VASP calculations.

The somewhat "primitive" way is creating a directory, writing up or copying from somewhere else the input files, most of the time INCAR, KPOINTS, POSCAR, POTCAR, and a job submission file depending on the supercomputer one is using, then submit the job and wait for the output files showing up in the same directory. 

But most of the time, one single run is far from enough to extract the information one wish to obtain. Rather, a series of runs with one or more varying parameters are needed, a post analysis should be done, and a cumulative but terse output in the form of a file or direct screen output should be provided to the user. Examples are __cut-off energy and k-points convergence tests__, __lattice constant (equilibrium volume) and bulk modulus test__, __elastic constant test__, and some __exotic tests__ that require a __gradual change__ somewhere in the __POSCAR__, such as a __rotation of several chosen ions__. No one should expect a user to have the patience to change the varying parameter of every single run in the series and submit them one by one, extract information from each run by hand, finally fit the curve, conclude the coefficient and preferably plot them in a saved figure. In order to tackle this problem, there are two automated ways.

The __first__ is wrapping up the job preparation, execution and analysis trilogy in an installable package, for a high level scripting language, such as Python. A famous example is [ASE] (https://wiki.fysik.dtu.dk/ase/). One needs to install the Python package and read the documentation on the location and functionality of their defined classes and methods, write up a Python script for a series of runs, oftentimes along with the post analysis part and feed it to the batch job submission script. It has a lot of advantages, especially if one deals with different kinds of _ab initio_ calculators and very routine work that all of them share, because it hides the basic file manipulation procedures, which are different for different calculators, and let users work with its defined Python classes and methods.

The __second__ is providing a collection of scripts, written in a shell scripting language, like Bash, and accompanying them with scripts written in a programming scripting language, like Python, in the same directory, with careful use of file names and input arguments, and let users add this single directory into their PATH environmental variable. The scripts will be used as executables. This project aims to strengthen such an approach. This method is less systematic than the former, but it also has some obvious advantages comparitively.

* It provides more flexibility in terms of file and directory manipulation, because Bash, as a shell language, deals with them more directly than Python, if one doesn't bother to use more Python packages in order to replace Bash. The directory that ASE runs in will not hold the full output of a series of runs, but only the last one, in an __equilibrium volume and bulk modulus__ subroutine. But with Bash one can conveniently create subdirectories to hold the input and output of each run in an automated way, thus keeping records of everything for future use.

* With all of the output files present, one can do far more checks and data collection with some knowledge of command line tools such as grep, sed, awk, find, etc. These tools are considered more universal than the layered classes and methods defined in ASE Python package. The executable scripts in a single directory are grouped more straigh-forwardly, showing the users only task related commands with simple keywords. As far as my own experience, VASP is already a delicate machine that requires enough memorization, there is really less motivation to learn a completely new Pythonic scheme. If there are more straight-forward ways to achieve the same goal, I would turn my head to them.

* If one writes up a Python script describing a series of runs with ASE, and submit the whole thing, if it reaches the walltime, or one or more of the runs crash, or the script is faulty somewhere, which are something that always happen, the script cannot be finished, and it would be difficult to pick up what has been terminated unexpectedly. The scheme is somewhat fragile. With the existence of multiple subdirectories, however, one can give each run in a series an independent job submission script, and submit them selectively. If one run reaches the walltime or simply fails, one can easily identify it, change the input parameter and submit it again, and other runs won't be affected. If there are mistakes in the Bash script, one will discover them before submission, because chances are that those mistakes are directly output on the screen.

* One can also freely separate the input file preparation, job submission and post-analysis, just as implemented in this project by __Prepare.sh__, __Fire.sh__, and __Display.sh__, or combine the first two for daily use when feeling secure as __Prep-fire.sh__. For a single run, one can just use __Prep-fast.sh__. The overall directory structure will be very neat. This opens a door widely to runs that depend on the output of previous runs. It's just several lines of commands away with all the convenient tools in the linux shell.

* To save more space, one can choose not to output or delete the WAVECAR, CHG, and CHGCAR.

Based on above reasons, which of course will be less persuasive for an experience ASE user who is already familiar with its Pythonic scheme, or if ASE improves its management of directories, I present the scripts that are used by me on a daily basis for dealing with subroutines described in the beginning of this section, and more, including the usual need for density of states plots, band structure plots, Bader charge transfer analysis, etc. The list will go on and can be mixed with other people's scripts conveniently.

Due to the low level of "abstraction", it is assumed that the user of this script collection is familiar with simple Bash syntax, command line tools such as cat, grep, sed, awk, find and some Python. Flexibility comes with some price, so one may need to look at the code and change some lines. But since one doesn't have to be forced to face a formidable Python package with intricate definitions, it will be much easier. I'll try to turn the hard coded parts to automated detection and input command line arguments, and make appropriate comments on the subroutines as I progress. This is still an ongoing project, but you are welcome to fork it and modify the scripts to your own need.

-------------------
    HAVEN'T FINALIZED YET BELOW

Scripts structure
---------------------

### For single-element runs
* Prepare.sh entest / lctest / c11 cubic / ... - to create a list of runs
* Fire.sh entest / lctest / c11 / ... - to submit a list of jobs
* Display.sh entest / lctest / c11 / ... - to do the analysis of a list of runs and determine physical quantities


* Prep-fire.sh - a wrapper on Prepare.sh and Fire.sh
* Prep-fast.sh - a fast way to create a folder with necessary input files


* Elastic.sh prep-fire / disp-solve / solve - to do a set of elastic constant determination runs, built on Prepare.sh, Fire.sh and Display.sh
* Dos-bs.sh scrun / dosrun / bsrun / plot-tdos / plot-ldos / plot-bs - to carry out self-consistent runs, density of states runs, band structure runs, to plot the total and local density of states, and to plot the band structure
* Bader.sh prerun /bader - to perform Bader charge transfer analysis

### For multi-element runs
* ElementRun.sh - an interactive script to execute the above subroutines for multiple elements

### Useful scripts
* Cellinfo.sh - to show the important info from OUTCAR
* Swap-file-name.sh - to exchange the names of two files
* Delfiles.sh - to delete the output files and start over for failed runs
* Re-relax.sh - to copy CONTCAR to POSCAR and start over relaxation

General subroutines
---------------------

1. Create a directory for a specific compound, say PtN. This is the top working directory
2. Go into the directory and create folder INPUT, putting INCAR, KPOINTS, POSCAR, concatenated POTCAR and qsub.parallel into it
3. Customize the input files to prepare the run, including reserving the strings that need to be replaced as @R@, @N@, etc.
4. Go back to the top working directory. Use Prepare.sh xxx yyy to perform a set of runs, or Prep-fast.sh xxx to simply create a folder xxx with proper substitutions
5. Use Fire.sh xxx to execute a set of qsub batch commands, or simply go into the folder prepared by Prep-fast.sh and qsub qsub.parallel
6. Using Display.sh xxx [M/P (lctest only)] to get the (step, energy) pairs written into xxx_output.txt and other instructive information. For lctest and elastic constant runs the fitting results and plot figures are also written
7. Rename the folder created during the process under the top working directory to whatever looks neat and tidy, i.e. labeling them with numbers

### To perform a complete set of elastic constant runs

1. Under the top working directory, execute Elastic.sh xxx prep-fire, i.e. cubic
2. Elastic.sh xxx disp-solve after every batched job is finished
3. Go check the elastic/elastic_output.txt for the results

### To perform runs for multiple elements

1. Create a directory, i.e. TMN
2. Create a folder INPUT, putting in necessary input files with necessary @R@ and @N@
2. Go into the directory and execute Element-run.sh. Do as suggested by the script
