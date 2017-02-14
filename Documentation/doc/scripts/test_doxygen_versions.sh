#!/bin/bash
if [ "$1" == '--help' ]; then 
  echo "Usage: $0 <doxygen_1> [doxygen_2] "
  echo "Compares the output of doxygen_1 and doxygen_2, "
  echo "where  doxygen_1 and doxygen_2 are valid paths to doxygen executables."
  echo "If doxygen_2 is not specified, the master branch of doxygen will be cloned, built and used as doxygen_2."
  exit 0
fi
#build reference
PATH_TO_1="$1"
PATH_TO_2="$2"
IS_ARG2=1
if [ -z $PATH_TO_2 ]; then
  IS_ARG2=0
fi

if [ -z $PATH_TO_1 ] || [ $(basename $PATH_TO_1) != "doxygen" ] || [ ! -e $PATH_TO_1 ]; then
  echo "Please specify a valid path to a doxygen executable."
  echo "$0 --help for more information."
  exit 0
fi
mkdir ./build_doc
cd ./build_doc
cmake -DCGAL_GENERATE_XML=ON -DDOXYGEN_EXECUTABLE="$PATH_TO_1" ../..  &> /dev/null
make -j7 doc  &> /dev/null
cd ../
bash compare_testsuites.sh $PWD/build_doc/doc_output
mv ./doc_data ./doc_ref

#download and build doxygen_master
if [ $IS_ARG2 == 0 ] || [ $(basename $PATH_TO_2) != "doxygen" ] || [ ! -e $PATH_TO_2 ]; then
  echo "No path to doxygen master were detected. Cloning..."
  mkdir doxygen_master
  git clone https://github.com/doxygen/doxygen.git doxygen_master  &> /dev/null
  cd doxygen_master
  mkdir build
  cd build
  cmake ..  &> /dev/null
  make  &> /dev/null
  cd ../..
  PATH_TO_2=../doxygen_master/build/bin/doxygen
  echo "done."
fi
#build doc with doxygen master
rm -rf ./build_doc
mkdir build_doc
cd ./build_doc
cmake -DCGAL_GENERATE_XML=ON -DDOXYGEN_EXECUTABLE="$PATH_TO_2" ../..  &> /dev/null
make -j7 doc  &> /dev/null
cd ../
bash ./compare_testsuites.sh $PWD/build_doc/doc_output $PWD/doc_ref
#clean-up
rm -rf ./build_doc
rm -rf ./doc_ref
rm -rf ./doc_data
if [ $IS_ARG2 == 0 ]; then
  rm -rf ./doxygen_master
fi
