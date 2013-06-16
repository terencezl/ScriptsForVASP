#!/bin/bash
if [ $1 == scrun ]; then
    Prep-fast.sh sc-dos-bs
    scaling_factor=$(grep "Equilibrium scaling factor is" lctest/lctest_output.txt | head -1 | awk '{print $5}')
    rwig1=$(cd elastic/*/0.00; Cellinfo.sh rwigs |awk '{print $4}')
    rwig2=$(cd elastic/*/0.00; Cellinfo.sh rwigs |awk '{print $5}')
    cd sc-dos-bs
    sed -i "s/@R@/$scaling_factor/g" POSCAR
    sed -i "s/LCHARG =/#LCHARG =/g" INCAR
    sed -i "s/#NPAR =/NPAR =/g" INCAR
    sed -i "s/#RWIGS =/RWIGS = $rwig1 $rwig2/g" INCAR
    qsub qsub.parallel
elif [ $1 == dosrun ]; then
    cd sc-dos-bs
    mkdir scrun
    mv * scrun 2> /dev/null
    (cd scrun; mv INCAR KPOINTS POTCAR POSCAR CHGCAR qsub.parallel ..)
    sed -i "s/#ICHARG =/ICHARG =/g" INCAR
    sed -i "s/#LCHARG =/LCHARG =/g" INCAR
    sed -i "s/#ISMEAR = 0/ISMEAR = -5/g" INCAR
    sed -i '3c Gamma' KPOINTS
    sed -i '4c 15 15 15' KPOINTS
    qsub qsub.parallel
elif [ $1 == plotdos ]; then
    cd sc-dos-bs
    pathname=(${PWD//\// })
    metal=$(echo ${pathname[$(( ${#pathname[@]} - 3 ))]})
    cryst_struct=$(echo ${pathname[$(( ${#pathname[@]} - 2 ))]})
    _Plot-ldos.py DOSCAR $metal $cryst_struct '' ''
fi
