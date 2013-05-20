#!/bin/bash
# Use it under the successfully relaxed folder.
# Onlly support 2 elements. Specify the number of atoms as
# Cellinfo.sh N1 N2

cat POTCAR |grep US |grep -v TITEL |head -1
cat POTCAR |grep ENMAX |head -1
cat POTCAR |grep RWIGS |head -1
cat POTCAR |grep US |grep -v TITEL |tail -1
cat POTCAR |grep ENMAX |tail -1
cat POTCAR |grep RWIGS |tail -1
cat OUTCAR |grep 'volume of cell' -A 7 |tail -8

if [ $1 ] && [ $2 ]; then
    Vpcell=$(cat OUTCAR |grep 'volume of cell' |tail -1| awk '{print $5;}')
    r1=$(cat POTCAR |grep RWIGS |head -1 |awk '{print $6;}')
    r2=$(cat POTCAR |grep RWIGS |tail -1 |awk '{print $6;}')
    R1=$(echo "scale=4; e(l($Vpcell*3/4/3.14159/($1+$2*($r2/$r1)^3))/3)" | bc -l)
    R2=$(echo "scale=4; $r2/$r1*$R1" | bc -l)
    echo "You should use $R1 $R2 as your RWIGS in INCAR to get 100% filling."
fi
