#!/bin/bash

# Intensive test: test a parameterization method with all models in data folder
# Usage: test_method.sh parameterization-method boundary-parameterization
# Example: test_method.sh floater circle 
# Example: test_method.sh lscm 2pts 

echo "                                     *                                     "
echo "                                    ***                                    "
echo "***************************************************************************"

./test.sh "$1" "$2" blech 2>&1
./test.sh "$1" "$2" cone 2>&1
# ./test.sh "$1" "$2" crater 2>&1
./test.sh "$1" "$2" cyl 2>&1
./test.sh "$1" "$2" david 2>&1
./test.sh "$1" "$2" fandisk-top 2>&1
./test.sh "$1" "$2" foot 2>&1
./test.sh "$1" "$2" hand 2>&1
./test.sh "$1" "$2" mannequin-devil 2>&1
./test.sh "$1" "$2" maxplanck 2>&1
./test.sh "$1" "$2" nefertiti 2>&1
./test.sh "$1" "$2" one_ring 2>&1
./test.sh "$1" "$2" oni 2>&1
./test.sh "$1" "$2" rect 2>&1
# ./test.sh "$1" "$2" rotor 2>&1
./test.sh "$1" "$2" S3-pet41-fixed 2>&1
./test.sh "$1" "$2" strate2 2>&1
./test.sh "$1" "$2" three_peaks 2>&1
./test.sh "$1" "$2" vdi_punch 2>&1
./test.sh "$1" "$2" venus-loop 2>&1
./test.sh "$1" "$2" vw1-remeshed 2>&1

echo "***************************************************************************"
echo "                                    ***                                    "
echo "                                     *                                     "
