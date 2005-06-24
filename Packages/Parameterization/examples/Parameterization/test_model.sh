#!/bin/bash

# Intensive test: test all surface parameterization methods with 1 model in data folder
# Usage: test_model.sh model_name_without_extension
# Example: test_model.sh rotor

echo "                    ***********************************"
date

./test.sh uniform square opennl "$1"
./test.sh floater circle taucs "$1"
./test.sh conformal circle taucs "$1"
./test.sh authalic square taucs "$1"
./test.sh lscm 2pts taucs "$1"

date
echo "                    ***********************************"

