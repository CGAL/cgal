#!/bin/bash

# Utility script to run the surface parameterization example without knowing the executable name
# Usage: test.sh parameterization-method boundary-parameterization solver output-format source-file-root
# Example: test.sh floater circle opennl obj sphere966 
# Example: test.sh lscm 2pts taucs eps sphere966 

# Find executable name (different on Windows and Unix)
[ -f ./release/polyhedron_ex_parameterization.exe ] && PARAM_APPLICATION="./release/polyhedron_ex_parameterization.exe"
[ -x ./polyhedron_ex_parameterization ] && PARAM_APPLICATION="./polyhedron_ex_parameterization"

# Find source file in data or data/extras folders
[ -f data/extras/"$5".off ] && SOURCE_FILE="data/extras/"$5".off"
[ -f data/"$5".off ] && SOURCE_FILE="data/"$5".off"

# Create test folder (if needed)
[ -d test ] || mkdir test

# echo on
set -x

$PARAM_APPLICATION -t "$1" -b "$2" -s "$3" -o "$4" $SOURCE_FILE test/test_"$5"_"$1"_"$2"."$4" 2>&1

