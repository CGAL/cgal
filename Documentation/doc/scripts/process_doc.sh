#!/bin/bash
if [ "$1" == '--help' ]; then
  echo "Usage: $0 <doxygen 1.8.4> <doxygen 1.8.13> [publish_dir]"
  echo "Compares the output of doxygen 1.8.13 and doxygen master to the one from doxygen 1.8.4, of this CGAL version, "
  echo "publish_dir is the path to the dir where the testsuite results are kept"
  echo "$0 must be called from doc/scripts"
  exit 0
fi

mkdir -p doc_1_8_4
mkdir -p doc_1_8_13
mkdir -p doc_master

PATH_TO_1_8_4="$1"
PATH_TO_1_8_13="$2"
PUBLISH_DIR="$3"

DOXYGEN_1=$($PATH_TO_1_8_4 --version)
DOXYGEN_2=$($PATH_TO_1_8_13 --version)


#######################################
## download and build doxygen_master ##
#######################################
echo "downloading and building master"
git clone https://github.com/doxygen/doxygen.git doxygen_master  1> /dev/null
cd doxygen_master
git pull  https://github.com/lrineau/doxygen.git 1> /dev/null
MASTER_DESCRIBE=$(git describe --tags)
mkdir -p build
cd build
cmake ..  1> /dev/null
make -j$NB_CORES 1> /dev/null
cd ../.. #scripts
PATH_TO_MASTER="$PWD/doxygen_master/build/bin/doxygen"
echo "done."

echo "comparing versions 1.8.4 and 1.8.13"
bash -$- test_doxygen_versions.sh $PATH_TO_1_8_4 $PATH_TO_1_8_13 $PWD/doc_1_8_4 $PWD/doc_1_8_13 $PUBLISH_DIR
mv diff.txt diff1.txt

echo "comparing versions 1.8.4 and master"
bash -$- test_doxygen_versions.sh $PATH_TO_1_8_4 $PATH_TO_MASTER $PWD/doc_1_8_4 $PWD/doc_master $PUBLISH_DIR
mv diff.txt diff2.txt

#update overview
CGAL_NAME=$(cat cgal_version)
python ${PWD}/testsuite.py --output-dir1 $PWD/doc_1_8_4/doc_output/ --output-dir2 $PWD/doc_1_8_13/doc_output/ --doc-log-dir1 $PWD/doc_1_8_4/doc_log/ \
  --doc-log-dir2 $PWD/doc_1_8_13/doc_log/ --doc-log-dir-master $PWD/doc_master/doc_log/ \
  --publish $PUBLISH_DIR --diff1 $PWD/diff1.txt --diff2 $PWD/diff2.txt --master-dir $PWD/doc_master/doc_output/ \
  --cgal-version "$CGAL_NAME" --do-copy-results --version-to-keep 10 --doxygen-version1 "$DOXYGEN_1" --doxygen-version2 "$DOXYGEN_2" --master-describe "$MASTER_DESCRIBE"

#clean-up
rm -rf ./doc_1_8_4 ./doc_1_8_13 ./doc_master ./doxygen_master
rm ./diff1.txt ./diff2.txt ./cgal_version

