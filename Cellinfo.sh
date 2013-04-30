#!/bin/bash
# Use it under the successfully relaxed folder.
# Onlly support 2 elements with the same number of atoms. Specify the number of atoms as
# Cellinfo.sh N 

Vpcell=$(cat OUTCAR |grep volume |tail -1|cut -c 24-)
r1=$(cat POTCAR |grep RWIGS |head -1| cut -c 33-42)
r2=$(cat POTCAR |grep RWIGS |tail -1| cut -c 33-42)

R1=$(echo "scale=4; e(l($Vpcell*3/4/$1/3.14159/(1+($r2/$r1)^3))/3)" | bc -l)
R2=$(echo "scale=4; $r2/$r1*$R1" | bc -l)

cat POTCAR |grep US |grep -v TITEL |head -1
cat POTCAR |grep ENMAX |head -1
cat POTCAR |grep RWIGS |head -1
cat POTCAR |grep US |grep -v TITEL |tail -1
cat POTCAR |grep ENMAX |tail -1
cat POTCAR |grep RWIGS |tail -1
cat OUTCAR |grep volume -A 7|tail -8
echo "You should use $R1 $R2 as your RWIGS in INCAR to get 100% filling."
