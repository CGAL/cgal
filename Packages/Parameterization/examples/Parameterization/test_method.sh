#!/bin/bash

# Intensive test: test a parameterization method with all models in data folder
# Usage: test_method.sh parameterization-method boundary-parameterization
# Example: test_method.sh floater circle 
# Example: test_method.sh lscm 2pts 

echo "                                     *                                     "
echo "                                    ***                                    "
echo "***************************************************************************"

# ./test.sh "$1" "$2" cube
./test.sh "$1" "$2" holes
./test.sh "$1" "$2" mannequin-devil
./test.sh "$1" "$2" mask_cone
./test.sh "$1" "$2" nefertiti
./test.sh "$1" "$2" rotor
./test.sh "$1" "$2" sphere966
./test.sh "$1" "$2" three_peaks


echo "***************************************************************************"
echo "                                    ***                                    "
echo "                                     *                                     "
