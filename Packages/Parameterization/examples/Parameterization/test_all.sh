#!/bin/bash

# Intensive test: test all surface parameterization methods with all models in data folder
date
./test_method.sh uniform square opennl
date
./test_method.sh floater circle opennl
date
./test_method.sh conformal circle opennl
date
./test_method.sh authalic square opennl
date
./test_method.sh lscm 2pts opennl
date