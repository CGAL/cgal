#!/bin/bash

# Usage: test.sh parameterization-method boundary-parameterization source-file-root
# Example: test.sh floater circle blech 
# Example: test.sh lscm 2pts blech 

set -x

./release/polyhedron_ex_parameterization.exe -t "$1" -b "$2" -o obj data/"$3".off > test/test_"$3"_"$1"_"$2".obj
