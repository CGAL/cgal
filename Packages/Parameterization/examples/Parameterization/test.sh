#!/bin/bash

# Utility script to run the surface parameterization example without knowing the executable name
# Usage: test.sh parameterization-method boundary-parameterization source-file-root
# Example: test.sh floater circle blech 
# Example: test.sh lscm 2pts blech 

# Find executable name (different on Windows and Unix)
[ -f ./release/polyhedron_ex_parameterization.exe ] && PARAM_APPLICATION="./release/polyhedron_ex_parameterization.exe"
[ -x ./polyhedron_ex_parameterization ] && PARAM_APPLICATION="./polyhedron_ex_parameterization"

# Create test folder (if needed)
mkdir test >/dev/null 2>/dev/null

# echo on
set -x

$PARAM_APPLICATION -t "$1" -b "$2" -o obj data/"$3".off test/test_"$3"_"$1"_"$2".obj 2>&1

