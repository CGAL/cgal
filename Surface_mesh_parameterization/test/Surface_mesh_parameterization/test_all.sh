#!/bin/bash

# Intensive test: test all surface parameterization methods with all models in data folder

# Find executable name (different on Windows and Unix)
[ -f ./release/extensive_parameterization_test.exe ] && PARAM_APPLICATION="./release/extensive_parameterization_test.exe"
[ -x ./extensive_parameterization_test ] && PARAM_APPLICATION="./extensive_parameterization_test"
if [ -z "$PARAM_APPLICATION" ]; then
    echo "Cannot find extensive_parameterization_test executable"
    exit 1
fi

# run test (echo on)
set -x
$PARAM_APPLICATION data/*.off 2>&1

