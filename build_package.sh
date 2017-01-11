#!/bin/bash

EXAMPLES="$1/examples/$1"
TEST="$1/test/$1"
DEMO="$1/demo/$1"
ROOT="$PWD"
if [ -d "$ROOT/$EXAMPLES" ]
then
  cd $ROOT/$EXAMPLES
  mkdir -p build
  cd build
  cmake -DCGAL_DIR="$ROOT/build" ..
  make
fi
if [ -d "$ROOT/$TEST" ]
then
  cd $ROOT/$TEST
  mkdir -p build
  cd build
  cmake -DCGAL_DIR="$ROOT/build" ..
  make
fi
#if [ "$1" != Polyhedron ]
#then
#  cd $ROOT/$DEMO
#  mkdir -p build
#  cd build
#  cmake -DCGAL_DIR="$ROOT/build" ..
#  make
#fi
