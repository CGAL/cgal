#!/bin/bash

set -e


if [ "$1" == '--help' ]; then
  echo "Usage: $0 <doxygen_1> <doxygen_2> <doc_1_dir> <doc_2_dir> <publish_dir>"
  echo "Compares the output of doxygen_1 and doxygen_2 of this CGAL version, "
  echo "where  doxygen_1 and doxygen_2 are valid paths to doxygen executables."
  echo "publish_dir is the path to the dir where the testsuite results are kept"
  echo "doc_##_dir is the path to the directory where the documentation will be output for the corresponding version."
  echo "$0 must be called from doc/scripts"
  exit 0
fi
#build reference
PATH_TO_1="$1"
PATH_TO_2="$2"
PUBLISH_DIR="$5"
BUILD_DIR_1="$3"
BUILD_DIR_2="$4"
NB_CORES="$(grep -c ^processor /proc/cpuinfo)"

if [ -z $PATH_TO_1 ] || [ ! -e $PATH_TO_1 ]; then
  echo "Please specify a valid path to a doxygen executable: $PATH_TO_1 is not valid."
  echo "$0 --help for more information."
  exit 0
fi

if [ ! -d $BUILD_DIR_1 ] || [ ! -d $BUILD_DIR_2 ] || [ ! -d $PUBLISH_DIR ]; then
  echo "doc_1_dir, doc_2_dir and publish_dir must be directories."
  echo "$0 --help for more information."
  exit 0
fi

HAS_REF=1
TEST=$(ls $BUILD_DIR_1)
if [ -z "$TEST" ]; then
  HAS_REF=0
fi
#Find the CGAL directory. If there is a directory called Documentation, this is a branch build.
#Else it is from a release.
TEMP=$PWD
IS_RELEASE=1
cd $PWD/../..
if [ "$(basename $PWD)" = 'Documentation' ]; then
  ROOT=$PWD/..
  IS_RELEASE=0
else
  ROOT=$PWD
fi
cd $TEMP #scripts
if [ "$HAS_REF" -ne "1" ]; then

  ################################################################
  ## Build a first time with Doxygen_1 and create the txt files ##
  ################################################################

  echo "Building reference documentation..."
  mkdir -p ./build_doc
  cd ./build_doc
  cmake -DCGAL_DOC_MATHJAX_LOCATION:STRING=../../MathJax -DCGAL_DOC_RELEASE=ON -DCGAL_GENERATE_XML=ON -DDOXYGEN_EXECUTABLE="$PATH_TO_1" ../..  1>> ./build_logs
  make -j$NB_CORES doc  &>> ./build_logs
  echo "done."
  cd ../ #scripts
  echo "Creating text files for diff...."
  bash -$- compare_testsuites.sh $PWD/build_doc/doc_output 1> /dev/null
  mv ./doc_data ./doc_ref
  cp -r doc_ref first_doc_ref
  echo "done."
  mv ./build_doc ./doc_dir
else
  echo "There is already a reference. Not re-building."
fi

##################################################################
## build doc with Doxygen_2, create the txt files and diff ##
##################################################################
echo "Building second documentation..."
mkdir -p build_doc
cd ./build_doc
cmake -DCGAL_DOC_MATHJAX_LOCATION:STRING=../../MathJax -DCGAL_DOC_RELEASE=ON -DCGAL_GENERATE_XML=ON -DDOXYGEN_EXECUTABLE="$PATH_TO_2" ../..  1>> ./build_logs
make -j$NB_CORES doc  &>> ./build_logs
echo "done."
cd ../ #scripts
DOXYGEN_1=$($PATH_TO_1 --version)
DOXYGEN_2=$($PATH_TO_2 --version)
echo "Comparing results..."
if [ "$HAS_REF" -eq "1" ]; then
  bash -$- ./compare_testsuites.sh $PWD/build_doc/doc_output $PWD/first_doc_ref 1> /dev/null
else
  bash -$- ./compare_testsuites.sh $PWD/build_doc/doc_output $PWD/doc_ref 1> /dev/null
fi
echo "done."
#add post-processing
cd ./build_doc
echo "Adding postprocessing..."
make -j$NB_CORES doc_with_postprocessing  &>> ./build_logs
echo "done."
cd .. #scripts
mv ./build_doc/* $BUILD_DIR_2
rm $BUILD_DIR_2/CMakeCache.txt


if [ "$HAS_REF" -ne "1" ]; then
  #######################################################################################################################
  ## rebuild docs to get the right doc_tags without GENERATE_XML because it ignores the EXCLUDE_SYMBOLS, ##
  ## which disrupts the logs                                                                                           ##
  #######################################################################################################################
  rm -rf ./doc_dir
  cd $BUILD_DIR_1
  cmake -DCGAL_DOC_MATHJAX_LOCATION:STRING=../../MathJax -DCGAL_DOC_RELEASE=ON -DCGAL_DOC_CREATE_LOGS="true" -DDOXYGEN_EXECUTABLE="$PATH_TO_1" ../..  1>> ./build_logs
  echo "Building reference documentation with postprocessing..."
  make -j$NB_CORES doc  &>> ./build_logs
  make -j$NB_CORES doc  &>> ./build_logs
  make -j$NB_CORES doc_with_postprocessing &>> ./build_logs
  echo "done."
  if [ $IS_RELEASE = 0 ]; then
    cd $ROOT
    mkdir -p ./build && cd ./build
    cmake -DWITH_CGAL_Core=false -DWITH_CGAL_ImageIO=false -DWITH_CGAL_Qt6=false .. 1>> ./build_logs
    CGAL_NAME="$(cat $PWD/VERSION)"
    cd $ROOT
    rm -rf ./build
    cd $ROOT/Documentation/doc/scripts
  else
    CGAL_NAME="$(cat $ROOT/VERSION)"
    cd $ROOT/doc/scripts
  fi
  echo "$CGAL_NAME">cgal_version
else
  echo "There is already a reference. Not re-building."
  rm -rf ./first_doc_ref
fi
  cd $BUILD_DIR_2
  cmake -DCGAL_DOC_MATHJAX_LOCATION:STRING=../../MathJax -DCGAL_DOC_RELEASE=ON -DCGAL_DOC_CREATE_LOGS="true" -DDOXYGEN_EXECUTABLE="$PATH_TO_2" ../..  1>> ./build_logs
  echo "Building reference documentation with postprocessing..."
  make -j$NB_CORES doc  &>> ./build_logs
  make -j$NB_CORES doc  &>> ./build_logs
  make -j$NB_CORES doc_with_postprocessing &>> ./build_logs
  echo "done."
  cd .. #scripts
  #get VERSION's content


echo "cleaning up"
#clean-up
rm -rf ./build_doc
rm -rf ./doc_ref
rm -rf ./doc_data
