#!/usr/bin/env bash

test_type="$1"
directory_name=bader
fname="$directory_name"_output.txt
shift 1

while getopts ":d:mf" opt; do
    case $opt in
    d)
        directory_name=$OPTARG
        ;;
    m)
        is_submit=true
        ;;
    f)
        is_override=true
        ;;
    \?)
        echo "Invalid option: -$OPTARG" >&2
        exit 1
        ;;
    :)
        echo "Option -$OPTARG requires an argument." >&2
        exit 1
        ;;
  esac
done

if [[ "$test_type" == test ]]; then
    if [[ -d "$directory_name" && $(ls -A "$directory_name") ]]; then
        echo -n "$directory_name/ contains files or sub-directories. "
        if [[ $is_override ]]; then
            echo "Overriding..."
        else
            echo "Escaping..."
            exit 1
        fi
    fi
    Prepare.sh "$directory_name" -fa qbader.serial
    cd "$directory_name"
    sed -i '/LAECHG/c LAECHG = .TRUE.' INCAR
#    echo -e 'NGXF = 250\nNGYF = 250\nNGZF = 250' >> INCAR
    sed -i '/NGXF/c NGXF = 250' INCAR
    sed -i '/NGYF/c NGYF = 250' INCAR
    sed -i '/NGZF/c NGZF = 250' INCAR
    sed -i '/LCHARG/c LCHARG = .TRUE.' INCAR
    sed -i '/NSW/c NSW = 0' INCAR
    [[ $is_submit ]] && qsub qsub.parallel

elif [[ "$test_type" == analysis ]]; then
    if [[ -d "$directory_name" && $(ls -A "$directory_name") ]]; then
        cd "$directory_name"
    else
        echo "$directory_name/ does not exist!"
        exit 1
    fi
    [[ $is_submit ]] && qsub qbader.serial

else
    echo "Specify what you are going to test!" >&2
    exit 1
fi

exit 0