#!/bin/bash

# Intensive test: test all surface parameterization methods with all models in data folder
./test_method.sh uniform square opennl  | tee test_all.log
./test_method.sh floater circle taucs   | tee -a test_all.log
./test_method.sh conformal circle taucs | tee -a test_all.log
./test_method.sh authalic square taucs  | tee -a test_all.log
./test_method.sh lscm 2pts opennl       | tee -a test_all.log
