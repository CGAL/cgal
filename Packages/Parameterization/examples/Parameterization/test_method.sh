#!/bin/bash

# Intensive test: test a parameterization method with all models in data folder
# Usage: test_method.sh parameterization-method boundary-parameterization solver
# Example: test_method.sh floater circle opennl
# Example: test_method.sh lscm 2pts taucs

echo "                                     *                                     "
echo "                                    ***                                    "
echo "***************************************************************************"

./test.sh "$1" "$2" "$3" cube
./test.sh "$1" "$2" "$3" holes
./test.sh "$1" "$2" "$3" mannequin-devil
./test.sh "$1" "$2" "$3" mask_cone
./test.sh "$1" "$2" "$3" nefertiti
./test.sh "$1" "$2" "$3" rotor
./test.sh "$1" "$2" "$3" sphere966
./test.sh "$1" "$2" "$3" three_peaks


echo "***************************************************************************"
echo "                                    ***                                    "
echo "                                     *                                     "
