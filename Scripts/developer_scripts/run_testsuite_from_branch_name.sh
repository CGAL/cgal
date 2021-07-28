#!/bin/bash

#To run: $1 = name of the user
#        $2 = name of the branch
#        $3 = base ref name (master, 5.1.x, 5.2.x, etc...)
#        $4 = number of the PR


if uname | grep -q -i cygwin; then
  #Is supposed to ignore \r as eol character.
  export SHELLOPTS
  set -o igncr
fi
source ~/.autofilterrc
(
cd $CGAL_ROOT
USER_REPO=$1
BRANCH_NAME=$2
BASE_NAME=$3
PR_NUMBER=$4

cd ${CGAL_ROOT}
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
TMP_LIST=$(git diff --name-only HEAD cgal/$BASE_NAME |cut -s -d/ -f1 |sort -u | xargs -I {} ls -d {}/package_info 2>/dev/null  |cut -d/ -f1 |egrep -v Installation||true)
LIST_OF_PKGS=""
for PKG in $(ls) ; do
  if [ -f $PKG/package_info/$PKG/dependencies ]; then
    if [ -n "$(comm -12  <(echo "$LIST_OF_PKGS"|sort) <(cat $PKG/package_info/$PKG/dependencies|sort))" ]; then
      LIST_OF_PKGS="$LIST_OF_PKGS $PKG"
    fi
  fi
done
if [ -f ${CGAL_ROOT}/list_test_packages ]; then rm ${CGAL_ROOT}/list_test_packages; fi
if [ "$LIST_OF_PKGS" != "" ]; then
  for f in $LIST_OF_PKGS
  do
    echo "echo \"$f\"" >> ${CGAL_ROOT}/list_test_packages
    echo "echo \"${f}_Examples\"" >> ${CGAL_ROOT}/list_test_packages
    echo "echo \"${f}_Demo\"" >> ${CGAL_ROOT}/list_test_packages
  done
fi
#create the release from the branch
echo " Create release..."
CGAL_VERSION="CGAL-$(sed -E 's/#define CGAL_VERSION (.*\..*)-dev/\1/' <(grep "#define CGAL_VERSION " Installation/include/CGAL/version.h))-${PR_NUMBER}"
echo "CGAL_VERSION = ${CGAL_VERSION}"> log
cmake -DGIT_REPO=${CGAL_GIT_DIR}/cgal -DDESTINATION=${CGAL_ROOT}/CGAL-TEST -DPUBLIC=OFF -DTESTSUITE=ON -DCGAL_VERSION=${CGAL_VERSION} -P ${CGAL_GIT_DIR}/cgal/Scripts/developer_scripts/cgal_create_release_with_cmake.cmake | tee log
echo "done."
DEST=$(sed -E 's/.*CGAL-TEST\/(.*)/\1/' log);

cd ${CGAL_ROOT}

if [ -L CGAL-I ]; then rm CGAL-I; fi
ln -s $PWD/CGAL-TEST/$DEST CGAL-I
echo "starting testsuite..."

./autotest_cgal -c 

echo "finished."
)>${CGAL_ROOT}/autotest.log2 2>&1 & 

echo "exit."
exit 0
