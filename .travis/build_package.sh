#!/bin/bash
set -e
EXAMPLES="$1/examples/$1"
TEST="$1/test/$1"
DEMO="$1/demo/$1"
ROOT="$PWD"
if [ -d "$ROOT/$EXAMPLES" ]
then
  cd $ROOT/$EXAMPLES
  mkdir -p build
  cd build
  cmake -DCGAL_DIR="$ROOT/build" -DCMAKE_CXX_FLAGS_RELEASE="-DNDEBUG" ..
  make -j2
fi
if [ -d "$ROOT/$TEST" ]
then
  cd $ROOT/$TEST
  mkdir -p build
  cd build
  cmake -DCGAL_DIR="$ROOT/build" -DCMAKE_CXX_FLAGS_RELEASE="-DNDEBUG" ..
  make -j2
fi
if [ "$1" == "CHECK" ]
then
  #parse current matrix and check that no package has been forgotten
  old_IFS=$IFS
  IFS=$'\n'
  COPY=0
  MATRIX=()
  for LINE in $(cat "$PWD/packages.txt")
  do
        MATRIX+="$LINE "
  done

  PACKAGES=()
  cd ..
  for f in *
  do
    if [ -d  "$f/examples/$f" ] || [ -d  "$f/test/$f" ] || [ -d  "$f/demo/$f" ]
        then
                PACKAGES+="$f "
        fi
  done


  DIFFERENCE=$(echo ${MATRIX[@]} ${PACKAGES[@]} | tr ' ' '\n' | sort | uniq -u)
  IFS=$old_IFS
  if [ "${DIFFERENCE[0]}" != "" ]
  then
        echo "Some packages are missing in the matrix."
        exit 1
  fi
  echo "Matrix is up to date."
fi
  exit 0

#if [ "$1" != Polyhedron ]
#then
#  cd $ROOT/$DEMO
#  mkdir -p build
#  cd build
#  cmake -DCGAL_DIR="$ROOT/build" ..
#  make
#fi
