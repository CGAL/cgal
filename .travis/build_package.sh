#!/bin/bash
set -e
IFS=$' '
ROOT="$PWD/../"
for ARG in $(echo "$@")
do
  if [ "$ARG" == "CHECK" ]
	then
    zsh -x $ROOT/Scripts/developer_scripts/test_merge_of_branch HEAD

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
	        echo "The matrix and the actual package list differ : ."
					echo ${DIFFERENCE[*]}
	        exit 1
	  fi
	  echo "Matrix is up to date."
    exit 0
	fi
	EXAMPLES="$ARG/examples/$ARG"
	TEST="$ARG/test/$ARG"
	DEMO="$ARG/demo/$ARG"
	if [ -d "$ROOT/$EXAMPLES" ]
	then
	  cd $ROOT/$EXAMPLES
	  mkdir -p build
	  cd build
	  cmake -DCGAL_DIR="$ROOT/build" -DCMAKE_CXX_FLAGS_RELEASE="-DCGAL_NDEBUG" ..
	  make -j2
	fi
	if [ -d "$ROOT/$TEST" ]
	then
	  cd $ROOT/$TEST
	  mkdir -p build
	  cd build
	  cmake -DCGAL_DIR="$ROOT/build" -DCMAKE_CXX_FLAGS_RELEASE="-DCGAL_NDEBUG" ..
  	make -j2
	fi
	
	#if [ "$ARG" != Polyhedron ]
	#then
	#  cd $ROOT/$DEMO
	#  mkdir -p build
	#  cd build
	#  cmake -DCGAL_DIR="$ROOT/build" ..
	#  make
	#fi
done

# Local Variables:
# tab-width: 2
# End:
