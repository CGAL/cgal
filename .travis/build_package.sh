#!/bin/bash
set -e
[ -n "$CGAL_DEBUG_TRAVIS" ] && set -x

CXX_FLAGS="-DCGAL_NDEBUG -ftemplate-backtrace-limit=0"

function mytime {
  /usr/bin/time -f "Spend time of %C: %E (real)" "$@"
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
    mytime sudo bash .travis/install_openmesh.sh
  fi


  if [ "$ARG" = "CHECK" ]
  then
    cd .travis
    mytime ./generate_travis.sh --check
    cd ..
    IFS=$old_IFS
    mytime zsh $ROOT/Scripts/developer_scripts/test_merge_of_branch HEAD
    #test dependencies
    cd $ROOT
    mytime bash Scripts/developer_scripts/cgal_check_dependencies.sh --check_headers /usr/bin/doxygen

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
    mytime cmake -DCMAKE_INSTALL_PREFIX=install/ -DCGAL_BUILD_THREE_DOC=TRUE ..
    mytime make install
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
    mytime cmake -DCMAKE_INSTALL_PREFIX=../../install -DCGAL_BUILD_THREE_DOC=TRUE ..
    exit 0
  fi

  if [ "$ARG" = "Installation" ]
  then
  mkdir build_dir
  cd build_dir
  cmake -DWITH_tests=ON -DBUILD_TESTING=ON ..
  ctest -j2 -L CGAL_cmake_testsuite --output-on-failure
  cd ..
  rm -rf ./build_dir
  #==-- configure all CGAL with -DWITH_examples=ON -DWITH_demos=ON -DWITH_tests=ON, and then launch CTest on a few labels. --==
  mkdir config_dir
  cd config_dir
  cmake -DWITH_examples=ON -DWITH_demos=ON -DWITH_tests=ON -DBUILD_TESTING=ON ..
  ctest -j2 -L AABB_tree --output-on-failure
  cd ..
  rm -rf ./config_dir
    exit 0
  fi

  IFS=$old_IFS

  if [ -n "$TRAVIS_PULL_REQUEST_BRANCH" ]; then
    DO_IGNORE=FALSE
    . $ROOT/.travis/test_package.sh "$ROOT" "$ARG"
    echo "DO_IGNORE is $DO_IGNORE"
    if [ "$DO_IGNORE" = "TRUE" ]; then
      continue
    fi
  fi
  IFS=$' '
  mkdir -p build-travis
  cd build-travis
  WITHDEMOS=ON
  if [ "$ARG" = "Polyhedron" ]; then
    WITHDEMOS=OFF
  fi
  EXTRA_CXX_FLAGS=
  case "$CC" in
    clang*)
      EXTRA_CXX_FLAGS="-Werror=inconsistent-missing-override"
      ;;
  esac


  mytime cmake -DCMAKE_CXX_FLAGS="${CXX_FLAGS} ${EXTRA_CXX_FLAGS}" -DCGAL_DONT_OVERRIDE_CMAKE_FLAGS:BOOL=ON -DBUILD_TESTING=ON -DWITH_tests=ON -DWITH_examples=ON -DWITH_demos=$WITHDEMOS ..
  mytime ctest -j2 -L $ARG'([_][A-Z]|$)'  -E execution___of__ --output-on-failure
done
IFS=$old_IFS
# Local Variables:
# tab-width: 2
# sh-basic-offset: 2
# End:
