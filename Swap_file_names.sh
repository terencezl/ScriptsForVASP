#!/usr/bin/env bash

mv "$1" swap-file-name
mv "$2" "$1"
mv swap-file-name "$2"
