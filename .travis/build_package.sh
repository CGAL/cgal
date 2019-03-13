#!/bin/bash
set -e
[ -n "$CGAL_DEBUG_TRAVIS" ] && set -x

CXX_FLAGS="-DCGAL_NDEBUG -ftemplate-backtrace-limit=0"

function time {
  /usr/bin/time -f "Spend time of %C: %E (real)" "$@"
}

function build_examples {
  mkdir -p build-travis
  cd build-travis
  time cmake -DCGAL_DIR="/usr/local/lib/cmake/CGAL" -DCMAKE_CXX_FLAGS="${CXX_FLAGS}" ..
  time make -j2
}

function build_tests {
  build_examples
}

function build_demo {
  mkdir -p build-travis
  cd build-travis
  EXTRA_CXX_FLAGS=
  case "$CC" in
    clang*)
      EXTRA_CXX_FLAGS="-Werror=inconsistent-missing-override"
      ;;
  esac
  time cmake -DCGAL_DIR="/usr/local/lib/cmake/CGAL" -DCGAL_DONT_OVERRIDE_CMAKE_FLAGS:BOOL=ON -DCMAKE_CXX_FLAGS="${CXX_FLAGS} ${EXTRA_CXX_FLAGS}" ..
  time make -j2
}
old_IFS=$IFS
IFS=$' '
ROOT="$PWD/.."
for ARG in $(echo "$@")
do
#skip package maintenance
  if [ "$ARG" = "Maintenance" ]; then
    continue
  fi
cd $ROOT
#install openmesh only if necessary
  if [ "$ARG" = "CHECK" ] || [ "$ARG" = BGL ] || [ "$ARG" = Convex_hull_3 ] ||\
     [ "$ARG" = Polygon_mesh_processing ] || [ "$ARG" = Property_map ] ||\
     [ "$ARG" = Surface_mesh_deformation ] || [ "$ARG" = Surface_mesh_shortest_path ] ||\
     [ "$ARG" = Surface_mesh_simplification ]; then
    time sudo bash .travis/install_openmesh.sh
  fi


  if [ "$ARG" = "CHECK" ]
  then
    cd .travis
    time ./generate_travis.sh --check
    cd ..
    IFS=$old_IFS
    time zsh $ROOT/Scripts/developer_scripts/test_merge_of_branch HEAD
    #test dependencies 
    cd $ROOT
    time bash Scripts/developer_scripts/cgal_check_dependencies.sh --check_headers /usr/bin/doxygen

    cd .travis
  	#parse current matrix and check that no package has been forgotten

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
	    if [ -d  "$f/package_info/$f" ]
	        then
	                PACKAGES+="$f "
	        fi
	  done	
	
	  DIFFERENCE=$(echo ${MATRIX[@]} ${PACKAGES[@]} | tr ' ' '\n' | sort | uniq -u)
	  IFS=$' '
	  if [ "${DIFFERENCE[0]}" != "" ]
	  then
	        echo "The matrix and the actual package list differ : ."
					echo ${DIFFERENCE[*]}
            echo "You should run generate_travis.sh."
	        exit 1
	  fi
	  echo "Matrix is up to date."
    #check if non standard cgal installation works
    cd $ROOT
    mkdir build_test
    cd build_test
    time cmake -DCMAKE_INSTALL_PREFIX=install/ ..
    time make install
    # test install with minimal downstream example
    mkdir installtest
    cd installtest
    touch main.cpp
    mkdir build
    echo 'project(Example)' >> CMakeLists.txt
    echo 'set(PROJECT_SRCS ${PROJECT_SOURCE_DIR}/main.cpp)' >> CMakeLists.txt
    echo 'find_package(CGAL REQUIRED)' >> CMakeLists.txt
    echo 'add_executable(${PROJECT_NAME} ${PROJECT_SRCS})' >> CMakeLists.txt
    echo 'target_link_libraries(${PROJECT_NAME} CGAL::CGAL)' >> CMakeLists.txt
    echo '#include "CGAL/remove_outliers.h"' >> main.cpp
    cd build
    time cmake -DCMAKE_INSTALL_PREFIX=../../install ..
    cd ..
    exit 0
  fi
  IFS=$old_IFS

  if [ -n "$TRAVIS_PULL_REQUEST" ] && [ "$ARG" != Polyhedron_demo ]; then
    DO_IGNORE=FALSE
    . $ROOT/.travis/test_package.sh "$ROOT" "$ARG"
    echo "DO_IGNORE is $DO_IGNORE"
    if [ "$DO_IGNORE" = "TRUE" ]; then
      continue
    fi
  fi
  IFS=$' '
  EXAMPLES="$ARG/examples/$ARG"
  TEST="$ARG/test/$ARG" 
  DEMOS=$ROOT/$ARG/demo/*

  if [ -d "$ROOT/$EXAMPLES" ]
  then
    cd $ROOT/$EXAMPLES
    if [ -f ./CMakeLists.txt ]; then
      build_examples
    else
      for dir in ./*
      do
        if [ -f $dir/CMakeLists.txt ]; then
          cd $ROOT/$EXAMPLES/$dir
          build_examples
        fi
      done
    fi
  elif [ "$ARG" != Polyhedron_demo ]; then
    echo "No example found for $ARG"
  fi

  if [ -d "$ROOT/$TEST" ]
  then
    cd $ROOT/$TEST
    if [ -f ./CMakeLists.txt ]; then
      build_tests
    else
      for dir in ./*
      do
        if [ -f $dir/CMakeLists.txt ]; then
          cd $ROOT/$TEST/$dir
          build_tests
        fi
      done
    fi
  elif [ "$ARG" != Polyhedron_demo ]; then
    echo "No test found for $ARG"
  fi
  #Packages like Periodic_3_triangulation_3 contain multiple demos
  for DEMO in $DEMOS; do
    DEMO=${DEMO#"$ROOT"}
    echo $DEMO
  	#If there is no demo subdir, try in GraphicsView
    if [ ! -d "$ROOT/$DEMO" ] || [ ! -f "$ROOT/$DEMO/CMakeLists.txt" ]; then
     DEMO="GraphicsView/demo/$ARG"
    fi
	  if [ "$ARG" != Polyhedron ] && [ -d "$ROOT/$DEMO" ]
  	then
      cd $ROOT/$DEMO
      build_demo
    elif [ "$ARG" != Polyhedron_demo ]; then
      echo "No demo found for $ARG"
	  fi
  done
  if [ "$ARG" = Polyhedron_demo ]; then
    DEMO=Polyhedron/demo/Polyhedron
    cd "$ROOT/$DEMO"
    build_demo
  fi
done
IFS=$old_IFS
# Local Variables:
# tab-width: 2
# sh-basic-offset: 2
# End:
