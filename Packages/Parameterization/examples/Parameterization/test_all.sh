#!/bin/bash

# Intensive test: test all surface parameterization methods with all models in data folder
date
./test_method.sh uniform square
date
./test_method.sh floater circle
date
./test_method.sh conformal circle
date
./test_method.sh authalic square
date
./test_method.sh lscm 2pts
date