#!/usr/bin/env bash

function qsub_replacer {
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
    sed -i "/#PBS -N/c #PBS -N $qname" $1
    sed -i "/^cd/c cd $PWD" $1
}

directory_name="$1"
shift 1

while getopts ":a:f" opt; do
    case $opt in
    a)
        additional_qsub_file=$OPTARG
        ;;
    f)
        is_override=true
        echo "-f triggered overriding existing directory."
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

if [[ -d "$directory_name" && $(ls -A "$directory_name") ]]; then
    echo -n "Directory contains files or sub-directories."
    if [[ $is_override ]]; then
        echo "Overriding..."
    else
        echo "Escaping..."
        exit 1
    fi
fi

mkdir -p "$directory_name" 2> /dev/null
cp INPUT/* "$directory_name"/
cd "$directory_name"
qsub_replacer qsub.parallel
if [[ -n "$additional_qsub_file" ]]; then qsub_replacer "$additional_qsub_file"; fi
