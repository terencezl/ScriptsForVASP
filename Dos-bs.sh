#!/bin/bash
if [ $1 == scrun ]; then
    Prep-fast.sh sc-dos-bs/scrun
#    rwig1=$(cd equi-relax; Cellinfo.sh rwigs |awk '{print $4}')
#    rwig2=$(cd equi-relax; Cellinfo.sh rwigs |awk '{print $5}')
    cd sc-dos-bs/scrun
    sed -i "/LWAVE/c #LWAVE = .FALSE." INCAR
    sed -i "/LCHARG/c #LCHARG = .FALSE." INCAR
    sed -i "/NPAR/c NPAR = 1" INCAR
#    sed -i "/RWIGS/c RWIGS = $rwig1 $rwig2" INCAR
    sed -i "/ISPIN/c ISPIN = 1" INCAR
    sed -i "/NSW/c NSW = 0" INCAR
    sed -i "/LORBIT/c LORBIT = 11" INCAR
    qsub qsub.parallel

elif [ $1 == dosrun ]; then
    Prep-fast.sh sc-dos-bs/dosrun
    cd sc-dos-bs/dosrun
    cp ../scrun/{INCAR,KPOINTS,CHGCAR,WAVECAR} .
    sed -i "/ICHARG/c ICHARG = 11" INCAR
    sed -i "/ISMEAR/c ISMEAR = -5" INCAR
    sed -i "/NEDOS/c NEDOS = 1501" INCAR
    sed -i '4c 21 21 21' KPOINTS
    echo "LMAXMIX = 4" >> INCAR
    qsub qsub.parallel

elif [ $1 == bsrun ]; then
    Prep-fast.sh sc-dos-bs/bsrun
    cd sc-dos-bs/bsrun
    cp ../scrun/{INCAR,KPOINTS,CHGCAR,WAVECAR} .
    echo "You must manually change the KPOINTS file before submitting job!"
#    qsub qsub.parallel

elif [ $1 == plot-ldos ]; then
    cd sc-dos-bs
    _Plot-ldos.py

elif [ $1 == plot-tdos ]; then
    cd sc-dos-bs
    _Plot-tdos.py

elif [ $1 == bs ]; then
    cd sc-dos-bs
    _Plot-bs.py
fi
