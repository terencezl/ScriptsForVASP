#!/bin/bash
if [ $1 == scrun ]; then
    Prep-fast.sh sc-dos-bs
    scaling_factor=$(grep "Equilibrium scaling factor is" lctest/lctest_output.txt | head -1 | awk '{print $5}')
    rwig1=$(cd elastic/*/0.00*; Cellinfo.sh rwigs |awk '{print $4}')
    rwig2=$(cd elastic/*/0.00*; Cellinfo.sh rwigs |awk '{print $5}')
    cd sc-dos-bs
    sed -i "s/@R@/$scaling_factor/g" POSCAR
    sed -i "/LWAVE/c #LWAVE = .FALSE." INCAR
    sed -i "/LCHARG/c #LCHARG = .FALSE." INCAR
    sed -i "/NPAR/c NPAR = 1" INCAR
    sed -i "/RWIGS/c RWIGS = $rwig1 $rwig2" INCAR
    sed -i "/ISPIN/c ISPIN = 1" INCAR
    sed -i "/NSW/c NSW = 0" INCAR
    echo "LMAXMIX = 4" >> INCAR
    qsub qsub.parallel
elif [ $1 == dosrun ]; then
    cd sc-dos-bs
    mkdir scrun
    mv * scrun 2> /dev/null
    (cd scrun; cp INCAR KPOINTS POTCAR POSCAR CHGCAR WAVECAR qsub.parallel ..)
    sed -i "/ICHARG/c ICHARG = 11" INCAR
    sed -i "/ISMEAR/c ISMEAR = -5" INCAR
    sed -i "/NEDOS/c NEDOS = 1501" INCAR
    sed -i '4c 21 21 21' KPOINTS
    qsub qsub.parallel
elif [ $1 == plot-ldos ]; then
    cd sc-dos-bs
    pathname=(${PWD//\// })
    metal=$(echo ${pathname[$(( ${#pathname[@]} - 3 ))]})
    cryst_struct=$(echo ${pathname[$(( ${#pathname[@]} - 2 ))]})
    _Plot-ldos.py DOSCAR $metal $cryst_struct '' ''
elif [ $1 == plot-tdos ]; then
    cd sc-dos-bs
    pathname=(${PWD//\// })
    metal=$(echo ${pathname[$(( ${#pathname[@]} - 3 ))]})
    cryst_struct=$(echo ${pathname[$(( ${#pathname[@]} - 2 ))]})
    _Plot-tdos.py DOSCAR $metal $cryst_struct
fi
