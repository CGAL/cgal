#!/bin/bash

#Will cd $1 and test package named $2
#to find out if it or one of its dependencies has changed in the current branch

DO_IGNORE=FALSE
cd $1
if [ ! -f "$2/package_info/$2/dependencies" ];then 
  echo "No dependencies found for $2"
  bash Scripts/developer_scripts/cgal_check_dependencies.sh --check_headers /usr/bin/doxygen
  
  
  exit 1
fi
LIST_OF_FILES=$(git diff --name-only origin/master... |cut -d/ -f1 |uniq |sort)
LIST_OF_DEPS=$(cat "$2/package_info/$2/dependencies")
echo "$LIST_OF_DEPS"
for flie in $LIST_OF_DEPS
do
  [[ $LIST_OF_FILES =~ (^|[[:space:]])$flie($|[[:space:]]) ]] && return
done
echo "Package ignored because none of its dependencies has been modified."
DO_IGNORE=TRUE

