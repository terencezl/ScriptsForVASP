#!/bin/bash
# Use it under a folder that has an OUTCAR
# Cellinfo.sh
# Cellinfo.sh rwigs [N1,N2] (Specify the number of atoms)
# Cellinfo.sh cubic [a_c11+2c12,a_c11-c12,a_c44]

Vpcell=$(cat OUTCAR |grep 'volume of cell' |tail -1| awk '{print $5;}')

if [ -z $1 ]; then
    grep -A 8 TITEL OUTCAR | grep -E 'TITEL | RWIGS | ENMAX'
    cat OUTCAR |grep 'volume of cell' -A 7 |tail -8

elif [ $1 == rwigs ]; then
    r_input=$(cat OUTCAR |grep RWIGS |grep wigner |awk '{print $6;}')
    r_input=[$(echo $r_input | sed 's/ /,/')]
    ElasticCoeff.py $1 "$2" $Vpcell "$r_input"

elif [[ $1 == cubic || $1 == tetragonal || $1 == orthorhombic || $1 == hexagonal || $1 == trigonal || $1 == monoclinic || $1 == triclinic ]]; then
    ElasticCoeff.py $1 "$2" $Vpcell
fi
