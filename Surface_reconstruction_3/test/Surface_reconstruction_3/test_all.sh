#!/bin/bash

#
# Test the Poisson Delaunay Reconstruction method with all off models in data folder
#

# Find executable name (different on Windows and Unix)
[ -f ./VC/debug/poisson_reconstruction_test.exe ] && PARAM_APPLICATION="./VC/debug/poisson_reconstruction_test.exe"
[ -f ./VC/release/poisson_reconstruction_test.exe ] && PARAM_APPLICATION="./VC/release/poisson_reconstruction_test.exe"
[ -f ./VC/x64/debug/poisson_reconstruction_test.exe ] && PARAM_APPLICATION="./VC/x64/debug/poisson_reconstruction_test.exe"
[ -f ./VC/x64/release/poisson_reconstruction_test.exe ] && PARAM_APPLICATION="./VC/x64/release/poisson_reconstruction_test.exe"
[ -x ./poisson_reconstruction_test ] && PARAM_APPLICATION="./poisson_reconstruction_test"
if [ -z "$PARAM_APPLICATION" ]; then
    echo "Cannot find poisson_reconstruction_test executable"
    exit 1
fi

# run test (echo on)
set -x
$PARAM_APPLICATION data/*.off 2>&1
set +x

#
# Test the Normal Estimation methods with all xyz models in data folder
#

# Find executable name (different on Windows and Unix)
[ -f ./VC/debug/normal_estimation_test.exe ] && PARAM_APPLICATION="./VC/debug/normal_estimation_test.exe"
[ -f ./VC/release/normal_estimation_test.exe ] && PARAM_APPLICATION="./VC/release/normal_estimation_test.exe"
[ -f ./VC/x64/debug/normal_estimation_test.exe ] && PARAM_APPLICATION="./VC/x64/debug/normal_estimation_test.exe"
[ -f ./VC/x64/release/normal_estimation_test.exe ] && PARAM_APPLICATION="./VC/x64/release/normal_estimation_test.exe"
[ -x ./normal_estimation_test ] && PARAM_APPLICATION="./normal_estimation_test"
if [ -z "$PARAM_APPLICATION" ]; then
    echo "Cannot find normal_estimation_test executable"
    exit 1
fi

# run test (echo on)
set -x
$PARAM_APPLICATION data/*.xyz 2>&1
set +x

