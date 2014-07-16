#!/bin/bash

# Intensive test: test all surface parameterization methods with 1 model in data folder
# Usage: test_model.sh source-file-root
# Example: test_model.sh rotor

echo ""
echo "                    ************************"
echo ""

./test.sh barycentric square opennl eps "$1"
echo "                                -"
./test.sh floater circle opennl obj "$1"
echo "                                -"
./test.sh conformal square opennl eps "$1"
echo "                                -"
# Skip authalic/opennl test which is very unstable
# ./test.sh authalic circle opennl eps "$1"
# echo "                                -"
./test.sh lscm 2pts opennl obj "$1"

echo ""
echo "                    ************************"
echo ""

