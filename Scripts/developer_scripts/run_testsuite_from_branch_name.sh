#!/bin/bash

#To run: $1 = name of the user
#        $2 = name of the branch
#        $3 = base ref name (master, 5.1.x, 5.2.x, etc...)
#        $4 = if given, update the VERSION_NUMBER


if uname | grep -q -i cygwin; then
  #Is supposed to ignore \r as eol character.
  export SHELLOPTS
  set -o igncr
fi
source "$(dirname $0)/.autofilterrc"
USER_REPO=$1
BRANCH_NAME=$2
BASE_NAME=$3
cd ${CGAL_ROOT}
scp mgimeno@cgal.geometryfactory.com:public_html/test_suite/VERSION_NUMBER .
cd ${CGAL_GIT_DIR}
if [ ! -d cgal ]; then
  git clone --depth 1 --no-single-branch https://github.com/CGAL/cgal.git 
  cd cgal
  git remote rename origin cgal
  cd ..
fi
cd cgal
git fetch --depth 1 cgal
git remote add $USER_REPO https://github.com/$USER_REPO/cgal.git
git fetch --depth 1 $USER_REPO
git checkout $BRANCH_NAME
git reset --hard $USER_REPO/$BRANCH_NAME
#setup the list_test_packages
LIST_OF_PKGS=$(git diff --name-only HEAD cgal/$BASE_NAME |cut -s -d/ -f1 |sort -u | xargs -I {} ls -d {}/package_info 2>/dev/null  |cut -d/ -f1 |egrep -v Installation||true)
if [ "$LIST_OF_PKGS" != "" ]; then
  if [ -f ${CGAL_ROOT}/list_test_packages ]; then rm ${CGAL_ROOT}/list_test_packages; fi
  for f in $LIST_OF_PKGS
  do
    echo "echo \"$f\"" >> ${CGAL_ROOT}/list_test_packages
    echo "echo \"${f}_Examples\"" >> ${CGAL_ROOT}/list_test_packages
    echo "echo \"${f}_Demo\"" >> ${CGAL_ROOT}/list_test_packages
  done
fi
(
#create the release from the branch
echo " Create release..."
cmake -DGIT_REPO=${CGAL_GIT_DIR}/cgal -DDESTINATION=${CGAL_ROOT}/CGAL-TEST -DPUBLIC=OFF -DTESTSUITE=ON -DCGAL_VERSION=0.0-Ic-$(cat ${CGAL_ROOT}/VERSION_NUMBER) -P ${CGAL_GIT_DIR}/cgal/Scripts/developer_scripts/cgal_create_release_with_cmake.cmake | tee log
echo "done."
DEST=$(sed -E 's/.*CGAL-TEST\/(.*)/\1/' log);

cd ${CGAL_ROOT}

if [ -L CGAL-I ]; then rm CGAL-I; fi
ln -s $PWD/CGAL-TEST/$DEST CGAL-I
echo "starting testsuite..."

./autotest_cgal -c 
)>${CGAL_ROOT}/autotest.log2 2>&1 & 
echo "finished."
if [ -n "$4" ]; then
  cd ${CGAL_ROOT}
  V=$(cat VERSION_NUMBER)
  V=$(($V+1))
  echo $V > VERSION_NUMBER
  scp VERSION_NUMBER mgimeno@cgal.geometryfactory.com:public_html/test_suite/VERSION_NUMBER
fi

echo "exit."
exit 0
