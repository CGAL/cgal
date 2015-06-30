#!/bin/sh

if [ ! -a CMakeLists.txt ] 
  then
    echo "CMakeLists.txt does not exist, generating..."
    cgal_create_CMakeLists -b program_options:timer -c Core
fi

python compileBenchmarks.py -f testModels.txt -d _modeldata -t benchmark_table -o benchmark_plot -s 6062699 -r 1 5,55,5
#python compileBenchmarks.py -k epeck -f simpleModels.txt -d _modeldata -t epeck_table -s 9894710 -r 1 -n 1 -q 1
#python compileBenchmarks.py -f hugeModels.txt -d _modeldata -t huge_table -s 5937524 -r 1 -n 1 -q 1
