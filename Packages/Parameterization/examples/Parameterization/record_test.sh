#!/bin/bash

# Utility script to run the surface parameterization example and record the STDERR output
# Usage: test.sh parameterization-method boundary-parameterization source-file-root 
# Example: record_test.sh floater circle blech 

# Find executable name (different on Windows and Unix)
[ -f ./release/polyhedron_ex_parameterization.exe ] && PARAM_APPLICATION="./release/polyhedron_ex_parameterization.exe"
[ -x ./polyhedron_ex_parameterization ] && PARAM_APPLICATION="./polyhedron_ex_parameterization"

# echo on
set -x

$PARAM_APPLICATION -t "$1" -b "$2" data/"$3".off 2>&1

