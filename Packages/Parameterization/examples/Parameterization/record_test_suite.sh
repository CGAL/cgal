#!/bin/bash

# Intensive test: run a set of tests of the surface parameterization example and record the STDERR output
# Usage: record_test_suite.sh

echo "***************************************************************************"
echo -e "\n"

echo "*** BUG 1 ***"
./record_test.sh authalic square venus-loop
echo -e "\n"

echo "*** BUG 2 ***"
./record_test.sh authalic square strate2
echo -e "\n"

echo "*** BUG 4 ***"
./record_test.sh uniform square one_ring
echo -e "\n"

#echo "*** BUG 3 ***"
#./record_test.sh lscm 2pts crater
#echo -e "\n"

echo "***************************************************************************"
echo -e "\n"
