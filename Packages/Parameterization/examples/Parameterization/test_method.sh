#!/bin/bash

# Intensive test: test a parameterization method with all models in data folder
# Usage: test_method.sh parameterization-method boundary-parameterization
# Example: test_method.sh floater circle 
# Example: test_method.sh lscm 2pts 

echo "                                     *                                     "
echo "                                    ***                                    "
echo "***************************************************************************"

./test.sh "$1" "$2" blech
./test.sh "$1" "$2" cone
# ./test.sh "$1" "$2" crater
./test.sh "$1" "$2" cyl
./test.sh "$1" "$2" david
./test.sh "$1" "$2" fandisk-top
./test.sh "$1" "$2" foot
./test.sh "$1" "$2" hand
./test.sh "$1" "$2" mannequin-devil
./test.sh "$1" "$2" maxplanck
./test.sh "$1" "$2" nefertiti
./test.sh "$1" "$2" one_ring
./test.sh "$1" "$2" oni
./test.sh "$1" "$2" rect
# ./test.sh "$1" "$2" rotor
./test.sh "$1" "$2" S3-pet41-fixed
./test.sh "$1" "$2" strate2
./test.sh "$1" "$2" three_peaks
./test.sh "$1" "$2" vdi_punch
./test.sh "$1" "$2" venus-loop
./test.sh "$1" "$2" vw1-remeshed

echo "***************************************************************************"
echo "                                    ***                                    "
echo "                                     *                                     "
