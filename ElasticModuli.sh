#!/bin/bash
# Use: In the working directory of elastic consts
# S. K. R. Patil, S. V. Khare, B. R. Tuttle, J. K. Bording, and S. Kodambaka, Mechanical stability of possible structures of PtN investigated using first-principles calculations, PHYSICAL REVIEW B 73, 104118 2006, DOI: 10.1103/PhysRevB.73.104118
# ElasticSolve.sh cubic/tetragonal/... original/alternative

cd elastic
if [ $1 == cubic ]; then
    dir_list="c11+2c12 c11-c12 c44"
    c11=$(grep "C11 is" $elastic_output.txt | head -1 | awk '{print $3}')
    c12=$(grep "C12 is" $elastic_output.txt | head -1 | awk '{print $3}')
    c44=$(grep "C44 is" $elastic_output.txt | head -1 | awk '{print $3}')
    econst=$c11" "$c12" "$c44
elif [ $1 == tetragonal ]; then echo
    dir_list="c11 c33 c44 5c11-4c12-2c13+c33 c11+c12-4c13+2c33 c11+c12-4c13+2c33+2c66"
    c11=$(grep "C11 is" $elastic_output.txt | head -1 | awk '{print $3}')
    c33=$(grep "C33 is" $elastic_output.txt | head -1 | awk '{print $3}')
    c44=$(grep "C44 is" $elastic_output.txt | head -1 | awk '{print $3}')
    c12=$(grep "C12 is" $elastic_output.txt | head -1 | awk '{print $3}')
    c13=$(grep "C13 is" $elastic_output.txt | head -1 | awk '{print $3}')
    c66=$(grep "C66 is" $elastic_output.txt | head -1 | awk '{print $3}')
    econst=$c11" "$c33" "$c44" "$c12" "$c13" "$c66
elif [ $1 == orthorhombic ]; then echo
    dir_list="c11 c22 c33 c44 c55 c66 4c11-4c12-4c13+c22+2c23+c33 c11-4c12+2c13+4c22-4c23+c33 c11+2c12-4c13+c22-4c23+4c33"
elif [ $1 == hexagonal ]; then echo
elif [ $1 == trigonal ]; then echo
elif [ $1 == monoclinic ]; then echo
elif [ $1 == triclinic ]; then echo
fi

Vpcell=$(grep 'volume of cell' $(cut -d' ' -f1 <<< $dir_list)/0.00/OUTCAR | tail -1 | awk '{print $5;}')
N_atoms=$(echo $(sed -n 6p $(cut -d' ' -f1 <<< $dir_list)/0.00/POSCAR))
N_atoms=[${N_atoms// /,}]
econst=[${econst// /,}]
echo '' >> elastic_output.txt
_ElasticModuli_solver.py $1 Vpcell "$N_atoms" "$econst"  | tee -a elastic_output.txt
