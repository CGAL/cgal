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

DO_COMPARE=TRUE
PATH_TO_SCRIPTS=${PWD}

#######################################
## download and build doxygen_master ##
 #######################################
echo "downloading and building master"
git clone https://github.com/doxygen/doxygen.git doxygen_master  1> /dev/null
if [ $? -ne 0 ]; then
  echo " clone of doxygen failed"
  DO_COMPARE=FALSE
else
  cd doxygen_master
  git pull  https://github.com/lrineau/doxygen.git 1> /dev/null
fi
if [ $? -ne 0 ] || [ "$DO_COMPARE" = "FALSE" ]; then
  echo " pull of doxygen failed"
  DO_COMPARE=FALSE
else
  MASTER_DESCRIBE=$(git describe --tags)
  mkdir -p build
  cd build
  cmake ..  1> /dev/null
fi
if [ $? -ne 0 ] || [ "$DO_COMPARE" = "FALSE" ]; then
  echo " cmake of doxygen failed"
  DO_COMPARE=FALSE
else
  make -j$NB_CORES 1> /dev/null
fi
cd $PATH_TO_SCRIPTS #scripts
PATH_TO_MASTER="$PWD/doxygen_master/build/bin/doxygen"
echo "done."

echo "comparing versions 1.8.4 and 1.8.13"
bash -$- test_doxygen_versions.sh $PATH_TO_1_8_4 $PATH_TO_1_8_13 $PWD/doc_1_8_4 $PWD/doc_1_8_13 $PUBLISH_DIR
mv diff.txt diff1.txt

echo "comparing versions 1.8.4 and master"
if [ "$DO_COMPARE" = "TRUE" ]; then
  bash -$- test_doxygen_versions.sh $PATH_TO_1_8_4 $PATH_TO_MASTER $PWD/doc_1_8_4 $PWD/doc_master $PUBLISH_DIR
fi
if [ $? -ne 0 ] || [ "$DO_COMPARE" = "FALSE" ]; then
  DO_COMPARE=FALSE
  echo " test_doxygen_versions with master failed"
  mv build_doc/build_logs doc_master/
else
  mv diff.txt diff2.txt
fi
#update overview
CGAL_NAME=$(cat cgal_version)
if [ "$DO_COMPARE" = "TRUE" ]; then
  python ${PWD}/testsuite.py --output-dir1 $PWD/doc_1_8_4/doc_output/ --output-dir2 $PWD/doc_1_8_13/doc_output/ --doc-log-dir1 $PWD/doc_1_8_4/doc_log/ \
    --doc-log-dir2 $PWD/doc_1_8_13/doc_log/ --doc-log-dir-master $PWD/doc_master/doc_log/ \
    --publish $PUBLISH_DIR --diff1 $PWD/diff1.txt --diff2 $PWD/diff2.txt --master-dir $PWD/doc_master/doc_output/ \
    --cgal-version "$CGAL_NAME" --do-copy-results --version-to-keep 10 --doxygen-version1 "$DOXYGEN_1" --doxygen-version2 "$DOXYGEN_2" --master-describe "$MASTER_DESCRIBE"
else
  echo "NO MASTER"
  python ${PWD}/testsuite.py --output-dir1 $PWD/doc_1_8_4/doc_output/ --output-dir2 $PWD/doc_1_8_13/doc_output/ --doc-log-dir1 $PWD/doc_1_8_4/doc_log/ \
    --doc-log-dir2 $PWD/doc_1_8_13/doc_log/ --doc-log-dir-master $PWD/doc_master/ \
    --publish $PUBLISH_DIR --diff1 $PWD/diff1.txt \
    --cgal-version "$CGAL_NAME" --do-copy-results --version-to-keep 10 --doxygen-version1 "$DOXYGEN_1" --doxygen-version2 "$DOXYGEN_2"
fi
#clean-up
rm -rf ./doc_1_8_4 ./doc_1_8_13 ./doc_master #./doxygen_master
rm ./diff1.txt ./cgal_version
if [ -f ./diff2.txt ]; then
  rm ./diff2.txt
fi

