#!/bin/bash
# Rerun.sh

cd $1 || exit 1
mkdir relax
mv * relax 2> /dev/null
cp relax/{INCAR,KPOINTS,POTCAR,CONTCAR,qsub.parallel} .
mv CONTCAR POSCAR
sed -i "s/NSW =/#NSW =/g" INCAR
sed -i "s/IBRION =/#IBRION =/g" INCAR
sed -i "s/EDIFFG =/#EDIFFG =/g" INCAR
sed -i "s/#LAECHG =/LAECHG =/g" INCAR
sed -i "s/#NG/NG/g" INCAR
qsub qsub.parallel
