#!/bin/bash

# Intensive test: test all surface parameterization methods with all models in data folder
./test_method.sh uniform square
./test_method.sh floater circle
./test_method.sh conformal circle
./test_method.sh authalic square
./test_method.sh lscm 2pts
