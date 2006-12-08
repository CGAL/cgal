#!/bin/bash

# Utility script to run the surface parameterization example without knowing the executable name
# Usage: test.sh parameterization-method border-parameterization solver output-format source-file-root
# Example: test.sh floater circle opennl obj sphere966
# Example: test.sh lscm 2pts taucs eps sphere966

# Find executable name (different on Windows and Unix)
[ -f ./release/polyhedron_ex_parameterization.exe ] && PARAM_APPLICATION="./release/polyhedron_ex_parameterization.exe"
[ -x ./polyhedron_ex_parameterization ] && PARAM_APPLICATION="./polyhedron_ex_parameterization"
if [[ ! -f $PARAM_APPLICATION ]] ; then
    echo "Cannot find polyhedron_ex_parameterization application"
    exit 1;
fi

# Find source file in data or data/extras folders
[ -f data/extras/"$5".off ] && SOURCE_FILE="data/extras/"$5".off"
[ -f data/"$5".off ] && SOURCE_FILE="data/"$5".off"

# Create test folder (if needed)
[ -d test ] || mkdir test

# Remove destination file (if needed)
DESTINATION_FILE="test/test_$5_$1_$2.$4"
[ -f $DESTINATION_FILE ] && rm -f $DESTINATION_FILE

# echo on
set -x

$PARAM_APPLICATION -t "$1" -b "$2" -s "$3" -o "$4" $SOURCE_FILE $DESTINATION_FILE 2>&1

