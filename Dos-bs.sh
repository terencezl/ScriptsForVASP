#!/bin/bash
if [ $1 == scrun ]; then
    Prep-fast.sh sc-dos-bs/scrun
#    rwig1=$(cd equi-relax; Cellinfo.sh rwigs |awk '{print $4}')
#    rwig2=$(cd equi-relax; Cellinfo.sh rwigs |awk '{print $5}')
    cd sc-dos-bs/scrun
#    sed -i "/LWAVE/c LWAVE = .TRUE." INCAR
    sed -i "/LCHARG/c LCHARG = .TRUE." INCAR
#    sed -i "/ISPIN/c ISPIN = 2" INCAR
    sed -i "/NSW/c NSW = 0" INCAR
    qsub qsub.parallel

elif [ $1 == dosrun ]; then
    Prep-fast.sh sc-dos-bs/dosrun
    cd sc-dos-bs/dosrun
    cp ../scrun/{INCAR,POSCAR,POTCAR,CHGCAR,WAVECAR} .
    sed -i "/ISMEAR/c ISMEAR = -5" INCAR
    sed -i "/NEDOS/c NEDOS = 2101" INCAR
    sed -i "/ICHARG/c ICHARG = 11" INCAR
    echo "LORBIT = 10" >> INCAR
    echo "LMAXMIX = 4" >> INCAR
    sed -i '4c 21 21 21' KPOINTS
    qsub qsub.parallel

elif [ $1 == bsrun ]; then
    Prep-fast.sh sc-dos-bs/bsrun
    cd sc-dos-bs/bsrun
    cp ../scrun/{INCAR,POSCAR,POTCAR,CHGCAR,WAVECAR} .
    sed -i "/ICHARG/c ICHARG = 11" INCAR
    echo "LORBIT = 10" >> INCAR
    echo "LMAXMIX = 4" >> INCAR
    if [ -f KPOINTS-bs ]; then
        mv KPOINTS-bs KPOINTS
    else
        echo "You must manually change the KPOINTS file before submitting job!"
        exit 1
    fi
    qsub qsub.parallel

elif [ $1 == plot-ldos ]; then
    cd sc-dos-bs/dosrun
    _Plot-ldos.py $2 $3

elif [ $1 == plot-tdos ]; then
    cd sc-dos-bs/dosrun
    _Plot-tdos.py $2

elif [ $1 == bs ]; then
    cd sc-dos-bs/bsrun
    _Plot-bs.py $2
fi
