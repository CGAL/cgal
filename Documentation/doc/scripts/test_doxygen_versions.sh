#!/bin/bash
if [ "$1" == '--help' ]; then 
  echo "Usage: $0 <doxygen_1> [doxygen_2] [CGAL_VERSION][publish_dir]"
  echo "Compares the output of doxygen_1 and doxygen_2 of this CGAL version, "
  echo "where  doxygen_1 and doxygen_2 are valid paths to doxygen executables."
  echo "If doxygen_2 is not specified, the master branch of doxygen will be cloned, built and used as doxygen_2."
  echo "CGAL_VERSION is the name that will appear on the results page. Default is CGAL-{version in cmake generated VERSION}"
  echo "publish_dir is the path to the dir where the testsuite results are kept"
  exit 0
fi
#build reference
PATH_TO_1="$1"
PATH_TO_2="$2"
CGAL_VERSION=$3
IS_ARG2=1
PUBLISH_DIR="$4"
if [ -z $PATH_TO_2 ]; then
  IS_ARG2=0
fi

if [ -z $PATH_TO_1 ] || [ $(basename $PATH_TO_1) != "doxygen" ] || [ ! -e $PATH_TO_1 ]; then
  echo "Please specify a valid path to a doxygen executable."
  echo "$0 --help for more information."
  exit 0
fi
#if no cgal version is specified, ask cmake"
if [ -n CGAL_VERSION ]; then
  mkdir ./shallow_build
  cd ./shallow_build 
  cmake ../../../.. &> /dev/null
  CGAL_VERSION="CGAL-$(cat VERSION)"
  cd ..
  rm -rf ./shallow_build
fi
mkdir ./build_doc
cd ./build_doc
cmake -DCGAL_GENERATE_XML=ON -DCGAL_DOC_CREATE_LOGS=true -DDOXYGEN_EXECUTABLE="$PATH_TO_1" ../..  &> /dev/null
make -j7 doc  &> /dev/null
make -j7 doc  &> /dev/null
cd ../ #scripts
bash compare_testsuites.sh $PWD/build_doc/doc_output
mv ./doc_data ./doc_ref

#download and build doxygen_master
if [ $IS_ARG2 == 0 ] || [ $(basename $PATH_TO_2) != "doxygen" ] || [ ! -e $PATH_TO_2 ]; then
  echo "No second path detected. Cloning Doxygen master branch..."
  git clone https://github.com/doxygen/doxygen.git doxygen_master  &> /dev/null
  cd doxygen_master
  mkdir build
  cd build
  cmake ..  &> /dev/null
  make -j7 &> /dev/null
  cd ../.. #scripts
  PATH_TO_2="$PWD/doxygen_master/build/bin/doxygen"
  echo "done."
fi
#build doc with doxygen master
mv ./build_doc ./doc_dir
mkdir build_doc
cd ./build_doc
cmake -DCGAL_GENERATE_XML=ON -DDOXYGEN_EXECUTABLE="$PATH_TO_2" ../..  &> /dev/null
make -j7 doc  &> /dev/null
cd ../ #scripts
DOXYGEN_1=$($PATH_TO_1 --version)
DOXYGEN_2=$($PATH_TO_2 --version)
bash ./compare_testsuites.sh $PWD/build_doc/doc_output $PWD/doc_ref $CGAL_VERSION $DOXYGEN_1 $DOXYGEN_2

#update overview
python ./testsuite.py --output-dir $PWD/doc_dir/doc_output/ --doc-log-dir $PWD/doc_dir/doc_log/ --publish $PUBLISH_DIR --diff $PWD/diff.txt

#clean-up
rm -rf ./build_doc
rm -rf ./doc_dir
rm -rf ./doc_ref
rm -rf ./doc_data
if [ $IS_ARG2 == 0 ]; then
  rm -rf ./doxygen_master
fi
