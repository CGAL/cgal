#!/bin/bash

# Intensive test: test all surface parameterization methods with 1 model in data folder
# Usage: test_model.sh source-file-root
# Example: test_model.sh rotor

echo "                    ***********************************"
date

./test.sh uniform square opennl obj "$1"
./test.sh floater circle taucs obj "$1"
./test.sh conformal circle taucs obj "$1"
./test.sh authalic square taucs obj "$1"
./test.sh lscm 2pts taucs obj "$1"

date
echo "                    ***********************************"

