#!/bin/bash
# PrepRerun.sh TEST_TYPE relax/static

fname="$1""_output.txt"
cd $1 || exit 1
dir_list=$(ls -F)

if [ $1 == relax ]; then
for n in $dir_list 
do
    if [[ "$n" == */ ]]
    then
        cd $n
        mkdir Temp
        mv INCAR KPOINTS POTCAR qsub.parallel CONTCAR Temp/
        rm * 2>/dev/null
        mv Temp/* .
        rm -r Temp
        mv CONTCAR POSCAR
        sed -i 's/POTIM = /#POTIM = /g' INCAR
    fi
done
echo -e "Second relaxation run after the first relaxation:\n" >> $fname

elif [ $1 == static ]; then
for n in $dir_list 
do
    if [[ "$n" == */ ]]
    then
        cd $n
        mkdir Temp
        mv INCAR KPOINTS POTCAR qsub.parallel CONTCAR Temp/
        rm * 2>/dev/null
        mv Temp/* .
        rm -r Temp
        mv CONTCAR POSCAR
        sed -i 's/ISIF = /#ISIF = /g' INCAR
        sed -i 's/NSW = /#NSW = /g' INCAR
        sed -i 's/IBRION = /#IBRION = /g' INCAR
    fi
done
echo -e "Static run after relaxation:\n" >> $fname

fi
