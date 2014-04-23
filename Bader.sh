#!/bin/bash

function qbader_replacer {
    qname_1=${PWD##*/}
    PWD_2=${PWD%/*}
    qname_2=${PWD_2##*/}
    PWD_3=${PWD_2%/*}
    qname_3=${PWD_3##*/}
    PWD_4=${PWD_3%/*}
    qname_4=${PWD_4##*/}
    qname="T$qname_4$qname_3$qname_2$qname_1"
    if [[ $(echo $qname | wc -c) > 17 ]]; then
        qname="T$qname_3$qname_2$qname_1"
    fi
    sed -i s/@N@/$qname/g qbader.serial
    sed -i s%@R@%$PWD%g qbader.serial
}

if [ $1 == prerun ]; then
    Fast-prep.sh bader
    cd bader || exit 1
    echo 'LAECHG = .TRUE.' >> INCAR
    echo -e 'NGXF = 250\nNGYF = 250\nNGZF = 250' >> INCAR
    sed -i '/LCHARG/c LCHARG = .TRUE.' INCAR
    sed -i '/NSW/c NSW = 0' INCAR
    qsub qsub.parallel
    cd ..

elif [ $1 == bader ]; then
    cd bader || exit 1
    cp ../INPUT/qbader.serial .
    qbader_replacer
    qsub qbader.serial
    cd ..

fi
