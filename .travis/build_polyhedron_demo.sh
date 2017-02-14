#!/bin/bash
set -e
cd ..
ROOT="$PWD"
cd ./Polyhedron/demo/Polyhedron
mkdir ./build
cd ./build
cmake -DCGAL_DIR="$ROOT/build" -DCMAKE_CXX_FLAGS_RELEASE="-DCGAL_NDEBUG" ..
make -j2
