#!/bin/bash
# Use: In the top working directory

if [ $1 == cubic ]; then
    dir_list="c11+2c12 c11-c12 c44"
#    dir_list="c11-c12 c44"
elif [ $1 == tetragonal ]; then echo
    dir_list="c11 c33 c44 5c11-4c12-2c13+c33 c11+c12-4c13+2c33 c11+c12-4c13+2c33+2c66"
elif [ $1 == orthorhombic ]; then echo
    dir_list="c11 c22 c33 c44 c55 c66 4c11-4c12-4c13+c22+2c23+c33 c11-4c12+2c13+4c22-4c23+c33 c11+2c12-4c13+c22-4c23+4c33"
elif [ $1 == hexagonal ]; then echo
elif [ $1 == trigonal ]; then echo
elif [ $1 == monoclinic ]; then echo
elif [ $1 == triclinic ]; then echo
fi

if [ $2 == prerun ]; then
    mkdir elastic-sc-dos 2> /dev/null
    cd elastic-sc-dos
    cp -r ../INPUT .
    scaling_factor=$(grep "Equilibrium scaling factor is" ../lctest/lctest_output.txt | head -1 | awk '{print $5}')
    sed -i "s/@R@/$scaling_factor/g" INPUT/POSCAR
    sed -i 's/#NSW /NSW /g' INPUT/INCAR
    sed -i 's/#IBRION /IBRION /g' INPUT/INCAR
    sed -i 's/#EDIFFG /EDIFFG /g' INPUT/INCAR
    sed -i "s/LCHARG =/#LCHARG =/g" INPUT/INCAR
    sed -i "s/#NPAR =/NPAR =/g" INPUT/INCAR

    for n in $dir_list
    do
        Prep-fire.sh $n $1
    done

elif [ $2 == scrun ]; then
    cd elastic-sc-dos || exit 1
    
    for n in $dir_list
    do
        for i in $(ls -F $n)
        do
            if [[ "$i" == */ ]]; then ratio_list=$ratio_list" "${i%/}; fi
        done
        for i in $ratio_list
        do
            rwig1=$(cd $n/$i; Cellinfo.sh rwigs | awk '{print $4}')
            rwig2=$(cd $n/$i; Cellinfo.sh rwigs | awk '{print $5}')
            cd $n/$i
	        mkdir prerun
            mv * prerun 2> /dev/null
            (cd prerun; cp INCAR KPOINTS POTCAR CONTCAR CHGCAR qsub.parallel ..)
            sed -i "s/#RWIGS =/RWIGS = $rwig1 $rwig2/g" INCAR
            sed -i 's/NSW /#NSW /g' INCAR
            sed -i 's/IBRION /#IBRION /g' INCAR
            sed -i 's/EDIFFG /#EDIFFG /g' INCAR
	        mv CONTCAR POSCAR
            cd ../..
        done
        unset ratio_list
        Fire.sh $n
    done
elif [ $2 == dosrun ]; then
    cd elastic-sc-dos || exit 1

    for n in $dir_list
    do
        for i in $(ls -F $n)
        do
            if [[ "$i" == */ ]]; then ratio_list=$ratio_list" "${i%/}; fi
        done
        for i in $ratio_list
        do
            cd $n/$i
            mkdir scrun
            mv * scrun 2> /dev/null
            (cd scrun; cp INCAR KPOINTS POTCAR POSCAR CHGCAR qsub.parallel ..; mv prerun ..)
            sed -i "s/#ICHARG =/ICHARG =/g" INCAR
            sed -i "s/#LCHARG =/LCHARG =/g" INCAR
            sed -i "s/#ISMEAR = 0/ISMEAR = -5/g" INCAR
            sed -i "3c Gamma" KPOINTS
            sed -i '4c 15 15 15' KPOINTS
            cd ../..
        done
        unset ratio_list
        Fire.sh $n
    done
    
elif [ $2 == plotdos ]; then
    cd elastic-sc-dos || exit 1
    mkdir dos-plots
    pathname=(${PWD//\// })
    metal=$(echo ${pathname[$(( ${#pathname[@]} - 3 ))]})
    cryst_struct=$(echo ${pathname[$(( ${#pathname[@]} - 2 ))]})
    
    for n in $dir_list
    do
        for i in $(ls -F $n)
        do
            if [[ "$i" == */ ]]; then ratio_list=$ratio_list" "${i%/}; fi
        done
        for i in $ratio_list
        do
            cd $n/$i
            _Plot-ldos.py DOSCAR $metal $cryst_struct $n $i
            cp LDOS-$metal"N"-$cryst_struct-$n-$i ../../dos-plots
            cd ../..
        done
        unset ratio_list
    done

# template for temporary use
elif [ $2 == fire ]; then
    cd elastic-sc-dos || exit 1
    for n in $dir_list
    do
        for i in $(ls -F $n)
        do
            if [[ "$i" == */ ]]; then ratio_list=$ratio_list" "${i%/}; fi
        done
        for i in $ratio_list
        do
            cd $n/$i
            sed -i "s/#ISMEAR = 0/ISMEAR = -5/g" INCAR
            sed -i "3c Gamma" KPOINTS
            sed -i '4c 15 15 15' KPOINTS
            cd ../..
        done
        unset ratio_list
        Fire.sh $n
    done
fi
