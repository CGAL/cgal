#!/bin/bash

# Intensive test: test all surface parameterization methods with all models in data folder

echo "***************************************************************"
date
echo "***************************************************************"

./test_model.sh cube
./test_model.sh holes
./test_model.sh mannequin-devil
./test_model.sh mask_cone
./test_model.sh nefertiti
./test_model.sh rotor
./test_model.sh sphere966
./test_model.sh three_peaks

echo "***************************************************************"
date
echo "***************************************************************"
