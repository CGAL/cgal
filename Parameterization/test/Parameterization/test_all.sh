#!/bin/bash

# Intensive test: test all surface parameterization methods with all models in data folder

echo "***************************************************************"
date
echo "***************************************************************"

# Find executable name (different on Windows and Unix)
[ -f ./release/extensive_parameterization_test.exe ] && PARAM_APPLICATION="./release/extensive_parameterization_test.exe"
[ -x ./extensive_parameterization_test ] && PARAM_APPLICATION="./extensive_parameterization_test"

# echo on
set -x

$PARAM_APPLICATION data/*.off 2>&1

# echo off
set +x

echo "***************************************************************"
date
echo "***************************************************************"
