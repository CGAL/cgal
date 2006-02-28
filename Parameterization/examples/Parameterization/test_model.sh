#!/bin/bash

# Intensive test: test all surface parameterization methods with 1 model in data folder
# Usage: test_model.sh source-file-root
# Example: test_model.sh rotor

echo "                    ************************"

./test.sh barycentric square opennl eps "$1"
echo "                                -"
./test.sh floater circle opennl obj "$1"
echo "                                -"
./test.sh conformal circle taucs obj "$1"
echo "                                -"
./test.sh authalic square taucs obj "$1"
echo "                                -"
./test.sh lscm 2pts taucs eps "$1"

echo "                    ************************"

