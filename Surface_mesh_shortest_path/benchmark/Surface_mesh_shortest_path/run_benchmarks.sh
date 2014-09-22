#!/bin/sh

if [ ! -a CMakeLists.txt ] 
  then
    echo "CMakeLists.txt does not exist, generating..."
    cgal_create_CMakeLists -b program_options:timer
fi

python compilePlots.py testModels.txt _modeldata benchmark_table benchmark_plot 6062699