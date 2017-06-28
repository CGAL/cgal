#!/bin/bash
if [ "$1" == '--help' ]; then
  echo "Usage: $0 <doxygen_1> [doxygen_2] [publish_dir]"
  echo "Compares the output of doxygen_1 and doxygen_2 of this CGAL version, "
  echo "where  doxygen_1 and doxygen_2 are valid paths to doxygen executables."
  echo "If doxygen_2 is not specified, the master branch of doxygen will be cloned, built and used as doxygen_2."
  echo "publish_dir is the path to the dir where the testsuite results are kept"
  echo "$0 must be called from doc/scripts"
  exit 0
fi
#build reference
PATH_TO_1="$1"
PATH_TO_2="$2"
PUBLISH_DIR="$3"
NB_CORES="$(grep -c ^processor /proc/cpuinfo)"

if [ -z $PATH_TO_1 ] || [ $(basename $PATH_TO_1) != "doxygen" ] || [ ! -e $PATH_TO_1 ]; then
  echo "Please specify a valid path to a doxygen executable."
  echo "$0 --help for more information."
  exit 0
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
################################################################
## Build a first time with Doxygen_1 and create the txt files ##
################################################################

echo "Building reference documentation..."
mkdir -p ./build_doc
cd ./build_doc
cmake -DCGAL_GENERATE_XML=ON -DDOXYGEN_EXECUTABLE="$PATH_TO_1" ../..  1> /dev/null
make -j$NB_CORES doc  &> /dev/null
echo "done."
cd ../ #scripts
echo "Creating text files for diff...."
bash compare_testsuites.sh $PWD/build_doc/doc_output 1> /dev/null
mv ./doc_data ./doc_ref
echo "done."

#######################################
## download and build doxygen_master ##
#######################################

if [ -z $PATH_TO_2 ] || [ $(basename $PATH_TO_2) != "doxygen" ] || [ ! -e $PATH_TO_2 ]; then
  echo "No second path detected. Cloning Doxygen master branch..."
  git clone https://github.com/doxygen/doxygen.git doxygen_master  1> /dev/null
  cd doxygen_master
  MASTER_DESCRIBE=$(git describe --tags)
  mkdir -p build
  cd build
  cmake ..  1> /dev/null
  make -j$NB_CORES 1> /dev/null
  cd ../.. #scripts
  PATH_TO_2="$PWD/doxygen_master/build/bin/doxygen"
  echo "done."
fi
mv ./build_doc ./doc_dir
##################################################################
## build doc with doxygen master, create the txt files and diff ##
##################################################################
echo "Building second documentation..."
mkdir -p build_doc
cd ./build_doc
cmake -DCGAL_GENERATE_XML=ON -DDOXYGEN_EXECUTABLE="$PATH_TO_2" ../..  1> /dev/null
make -j$NB_CORES doc  &> /dev/null
echo "done."
cd ../ #scripts
DOXYGEN_1=$($PATH_TO_1 --version)
DOXYGEN_2=$($PATH_TO_2 --version)
echo "Comparing results..."
bash ./compare_testsuites.sh $PWD/build_doc/doc_output $PWD/doc_ref 1> /dev/null
echo "done."
#add post-processing
cd ./build_doc
echo "Adding postprocessing..."
make -j$NB_CORES doc_with_postprocessing  &> /dev/null
echo "done."
cd .. #scripts
mv ./build_doc ./doc_master


#######################################################################################################################
## rebuild doc with Doxygen_1 to get the right doc_tags without GENERATE_XML because it ignores the EXCLUDE_SYMBOLS, ##
## which disrupts the logs                                                                                           ##
#######################################################################################################################
rm -rf ./doc_dir
mkdir ./doc_dir
cd ./doc_dir
cmake -DCGAL_DOC_CREATE_LOGS="true" -DDOXYGEN_EXECUTABLE="$PATH_TO_1" ../..  1> /dev/null
echo "Building reference documentation with postprocessing..."
make -j$NB_CORES doc  &> /dev/null
make -j$NB_CORES doc  &> /dev/null
make -j$NB_CORES doc_with_postprocessing &> /dev/null
echo "done."
cd .. #scripts
#get VERSION's content
if [ $IS_RELEASE = 0 ]; then
  cd $ROOT
  mkdir -p ./build && cd ./build
  cmake .. 1> /dev/null
  CGAL_NAME="$(cat $PWD/VERSION)"
  cd $ROOT
  rm -rf ./build
  cd $ROOT/Documentation/doc/scripts
else
  CGAL_NAME="$(cat $ROOT/VERSION)"
fi

#update overview

python ./testsuite.py --output-dir $PWD/doc_dir/doc_output/ --doc-log-dir $PWD/doc_dir/doc_log/ \
  --publish $PUBLISH_DIR --diff $PWD/diff.txt --master-dir $PWD/doc_master/doc_output/ \
  --cgal-version "$CGAL_NAME" --do-copy-results --version-to-keep 10 --doxygen-version "$DOXYGEN_1" --master-describe "$MASTER_DESCRIBE"

#clean-up
rm -rf ./build_doc
rm -rf ./doc_dir
rm -rf ./doc_master
rm -rf ./doc_ref
rm -rf ./doc_data
if [ -z $PATH_TO_2 ]; then
  rm -rf ./doxygen_master
fi
