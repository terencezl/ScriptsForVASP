#!/usr/bin/env bash
# General view and calculation of cell parameters, including rwigs
# Use it under a folder that has an OUTCAR and that has the equilibrium state of ions
# CellInfo.sh
# CellInfo.sh rwigs

test_type="$1"
Vpcell=$(cat OUTCAR |grep 'volume of cell' |tail -1| awk '{print $5;}')

if [[ -z "$test_type" ]]; then
    grep -A 8 TITEL OUTCAR | grep -E 'TITEL | RWIGS | ENMAX'
    cat OUTCAR |grep 'volume of cell' -A 7 |tail -8

elif [[ "$test_type" == rwigs ]]; then
    N_atoms=$(echo $(sed -n 6p POSCAR))
    N_atoms=[${N_atoms// /,}]
    r_input=$(echo $(cat OUTCAR | grep RWIGS | grep wigner | awk '{print $6;}'))
    r_input=[${r_input// /,}]
    _cell_info_solver.py $1 $Vpcell "$N_atoms" "$r_input"

else
    echo "Specify what you are going to show!" >&2
    exit 1
fi
