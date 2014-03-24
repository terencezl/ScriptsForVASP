#!/bin/bash

cd "$1" || exit 1
mkdir Temp
mv INCAR KPOINTS POTCAR POSCAR WAVECAR CHGCAR CONTCAR qsub.parallel Temp/
rm * 2>/dev/null
mv Temp/* .
#mv CONTCAR POSCAR
rm -r Temp
