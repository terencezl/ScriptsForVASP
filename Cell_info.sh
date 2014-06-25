#!/bin/bash
# General view and calculation of cell parameters, including rwigs
# Use it under a folder that has an OUTCAR and that has the equilibrium state of ions
# Cell_info.sh
# Cell_info.sh rwigs

Vpcell=$(cat OUTCAR |grep 'volume of cell' |tail -1| awk '{print $5;}')

if [ -z $1 ]; then
    grep -A 8 TITEL OUTCAR | grep -E 'TITEL | RWIGS | ENMAX'
    cat OUTCAR |grep 'volume of cell' -A 7 |tail -8

elif [ $1 == rwigs ]; then
    N_atoms=$(echo $(sed -n 6p POSCAR))
    N_atoms=[${N_atoms// /,}]
    r_input=$(echo $(cat OUTCAR | grep RWIGS | grep wigner | awk '{print $6;}'))
    r_input=[${r_input// /,}]
    _cell_info_solver.py $1 $Vpcell "$N_atoms" "$r_input"
fi
