#!/bin/bash

cd electronic 2>/dev/null
if [ $1 == scrun ]; then
    mkdir electronic && cd electronic
    cp -r ../INPUT .
    sed -i "/ISMEAR/c ISMEAR = 1" INPUT/INCAR
    sed -i "/NSW/c NSW = 0" INPUT/INCAR
    sed -i "/LCHARG/c LCHARG = .TRUE." INPUT/INCAR
    sed -i "/LMAXMIX/c LMAXMIX = 4" INPUT/INCAR
    Fast-prep.sh scrun
    cd scrun
    qsub qsub.parallel

elif [ $1 == dosrun ]; then
    Fast-prep.sh dosrun
    cd dosrun
    cp -l ../scrun/CHGCAR .
    cp ../scrun/INCAR .
    sed -i "/ISMEAR/c ISMEAR = -5" INCAR
    sed -i "/NEDOS/c NEDOS = 2101" INCAR
    sed -i "/ICHARG/c ICHARG = 11" INCAR
    if [[ $2 == rwigs ]]; then
        rwigs=$(cd ../scrun; Cellinfo.sh rwigs |awk '{print $4}')
        sed -i "/RWIGS/c RWIGS = ${rwigs//,/ }" INCAR
        sed -i "/NPAR/c NPAR = 1" INCAR
        sed -i "/LORBIT/c LORBIT = 5"  INCAR
    else
        sed -i "/LORBIT/c LORBIT = 11"  INCAR
    fi

    sed -i '4c 21 21 21' KPOINTS
    qsub qsub.parallel

elif [ $1 == bsrun ]; then
    Fast-prep.sh bsrun
    cd bsrun
    cp -l ../scrun/CHGCAR .
    cp ../scrun/INCAR .
    sed -i "/ICHARG/c ICHARG = 11" INCAR
    sed -i "/LORBIT/c LORBIT = 11" INCAR
    if [ -f KPOINTS-bs ]; then
        mv KPOINTS-bs KPOINTS
    else
        echo "You must manually change the KPOINTS file before submitting job!"
        exit 1
    fi
    qsub qsub.parallel

elif [ $1 == lobster-kp ]; then
    Fast-prep.sh lobster-kp
    cd lobster-kp
    sed -i "/ISYM/c ISYM = 0" INCAR
    sed -i "/LSORBIT/c LSORBIT = .TRUE." INCAR
    sed -i "/ISMEAR/c ISMEAR = -5" INCAR
    sed -i '4c 17 17 17' KPOINTS
    qsub qsub.parallel.lobster

elif [ $1 == lobster ]; then
    Fast-prep.sh lobster
    cd lobster
    cp -l ../scrun/CHGCAR .
    cp ../scrun/INCAR .
    cp ../lobster-kp/IBZKPT KPOINTS
    sed -i "/ISMEAR/c ISMEAR = -5" INCAR
    sed -i "/NBANDS/c NBANDS = 60" INCAR
    sed -i "/NEDOS/c NEDOS = 2101" INCAR
    sed -i "/ICHARG/c ICHARG = 11" INCAR
    sed -i "/LORBIT/c LORBIT = 11"  INCAR
    qsub qsub.parallel

elif [ $1 == plot-ldos ]; then
    cd dosrun
    _Plot-ldos.py $2 $3

elif [ $1 == plot-tdos ]; then
    cd dosrun
    _Plot-tdos.py $2

elif [ $1 == bs ]; then
    cd bsrun
    _Plot-bs.py $2
fi
