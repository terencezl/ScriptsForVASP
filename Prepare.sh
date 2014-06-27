#!/usr/bin/env bash

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
_qsub_replacer.sh qsub.parallel
if [[ -n "$additional_qsub_file" ]]; then _qsub_replacer.sh "$additional_qsub_file"; fi
