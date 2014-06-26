#!/usr/bin/env bash

cd "$1" || exit 1
#mkdir Temp
#mv INCAR KPOINTS POTCAR POSCAR WAVECAR CHGCAR CONTCAR qsub.parallel Temp/
#rm * 2>/dev/null
#mv Temp/* .
#rm -r Temp

rm XDATCAR PCDAT CONTCAR CHG WAVECAR EIGENVAL PROCAR vasprun.xml OUTCAR OSZICAR DOSCAR *.o*
