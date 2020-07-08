#!/bin/bash
set -e
[ -n "$CGAL_DEBUG_TRAVIS" ] && set -x

CXX_FLAGS="-DCGAL_NDEBUG -ftemplate-backtrace-limit=0"

function mytime {
  /usr/bin/time -f "Spend time of %C: %E (real)" "$@"
}

function build_examples {
  mkdir -p build-travis
  cd build-travis
  mytime cmake -DCGAL_DIR="/usr/local/lib/cmake/CGAL" -DCMAKE_CXX_FLAGS="${CXX_FLAGS}" -DCGAL_BUILD_THREE_DOC=TRUE ..
  mytime make -j2 VERBOSE=1
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
  mytime cmake -DCGAL_DIR="/usr/local/lib/cmake/CGAL" -DCGAL_DONT_OVERRIDE_CMAKE_FLAGS:BOOL=ON -DCMAKE_CXX_FLAGS="${CXX_FLAGS} ${EXTRA_CXX_FLAGS}"  ..
  mytime make -j2 VERBOSE=1
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
    exit 0
  fi

  if [ "$ARG" = "Installation" ]
  then

# ==-- CGAL installed in a non-default prefix, with the CGAL_Qt5 header-files. --==
    cd $ROOT
    mkdir build_test
    cd build_test
    mytime cmake -DCMAKE_INSTALL_PREFIX=install/ ..
    mytime make install
    # test install with minimal downstream example
    mkdir installtest
    cd installtest
    touch main.cpp
    mkdir build
    # https://doc.cgal.org/latest/Triangulation_2/Triangulation_2_2draw_triangulation_2_8cpp-example.html and  https://doc.cgal.org/latest/Algebraic_foundations/Algebraic_foundations_2interoperable_8cpp-example.html
    cat << EOF > main.cpp
      #include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
      #include <CGAL/Triangulation_2.h>
      #include <CGAL/draw_triangulation_2.h>
      #include <CGAL/basic.h>
      #include <CGAL/Coercion_traits.h>
      #include <CGAL/IO/io.h>
      #include <fstream>
      typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
      typedef CGAL::Triangulation_2<K>                            Triangulation;
      typedef Triangulation::Point                                Point;

      template <typename A, typename B>
      typename CGAL::Coercion_traits<A,B>::Type
      binary_func(const A& a , const B& b){
          typedef CGAL::Coercion_traits<A,B> CT;
          CGAL_static_assertion((CT::Are_explicit_interoperable::value));
          typename CT::Cast cast;
          return cast(a)*cast(b);
      }

      int main(int argc, char**) {
        std::cout<< binary_func(double(3), int(5)) << std::endl;
        std::cout<< binary_func(int(3), double(5)) << std::endl;
        std::ifstream in("data/triangulation_prog1.cin");
        std::istream_iterator<Point> begin(in);
        std::istream_iterator<Point> end;
        Triangulation t;
        t.insert(begin, end);
        if(argc == 3) // do not test Qt5 at runtime
          CGAL::draw(t);
        return EXIT_SUCCESS;
       }
EOF
    cat << EOF > "CMakeLists.txt"
      cmake_minimum_required(VERSION 3.1...3.15)
      find_package(CGAL COMPONENTS Qt5)
      add_definitions(-DCGAL_USE_BASIC_VIEWER -DQT_NO_KEYWORDS)
      add_executable(test main.cpp)
      target_link_libraries(test PUBLIC CGAL::CGAL_Qt5)
EOF

    cd build
    mytime cmake -DCMAKE_INSTALL_PREFIX=../../install ..
    make -j2
    ./test
    cd ..
    rm -rf ./build
    cd ..
    rm -rf ./install/*

# ==-- CGAL installed in a non-default prefix, without the CGAL_Qt5 header-files. --==
     # https://doc.cgal.org/latest/Triangulation_2/Triangulation_2_2draw_triangulation_2_8cpp-example.html and  https://doc.cgal.org/latest/Algebraic_foundations/Algebraic_foundations_2interoperable_8cpp-example.html
    mytime cmake -DCMAKE_INSTALL_PREFIX=install/ -DCGAL_BUILD_THREE_DOC=TRUE -DWITH_CGAL_Qt5=OFF ..
    mytime make install
    # test install with minimal downstream example
    cd installtest
    touch main.cpp
    mkdir build
    cat << EOF > "CMakeLists.txt"
      cmake_minimum_required(VERSION 3.1...3.15)
      find_package(CGAL)
      add_definitions(-DQT_NO_KEYWORDS)
      add_executable(test main.cpp)
      target_link_libraries(test PUBLIC CGAL::CGAL)
EOF
    cat << EOF > main.cpp
    #include <iostream>
    #include <CGAL/Simple_cartesian.h>
    typedef CGAL::Simple_cartesian<double> Kernel;
    typedef Kernel::Point_2 Point_2;
    typedef Kernel::Segment_2 Segment_2;
    int main()
    {
      Point_2 p(1,1), q(10,10);
      std::cout << "p = " << p << std::endl;
      std::cout << "q = " << q.x() << " " << q.y() << std::endl;
      std::cout << "sqdist(p,q) = "
                << CGAL::squared_distance(p,q) << std::endl;

      Segment_2 s(p,q);
      Point_2 m(5, 9);

      std::cout << "m = " << m << std::endl;
      std::cout << "sqdist(Segment_2(p,q), m) = "
                << CGAL::squared_distance(s,m) << std::endl;
      std::cout << "p, q, and m ";
      switch (CGAL::orientation(p,q,m)){
      case CGAL::COLLINEAR:
        std::cout << "are collinear\n";
        break;
      case CGAL::LEFT_TURN:
        std::cout << "make a left turn\n";
        break;
      case CGAL::RIGHT_TURN:
        std::cout << "make a right turn\n";
        break;
      }
      std::cout << " midpoint(p,q) = " << CGAL::midpoint(p,q) << std::endl;
      return 0;
    }
EOF
    cd build
    mytime cmake -DCMAKE_INSTALL_PREFIX=../../install ..
    make -j2
    ./test
    cd ..
    rm -rf ./build
    cd ../..
    rm -rf ./install_test

#==-- configure CGAL (as header-only), then use its build-directory as CGAL_DIR (we might want to support that use-case to be compatible with our usage with non-header-only, for CGAL-4.14) --==
mkdir config_dir
cd config_dir
mytime cmake  ..
cd ..
mkdir test_dir
#==-- test install with minimal downstream example --==
cd test_dir
touch main.cpp
mkdir build
cat << EOF > main.cpp
  #include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
  #include <CGAL/Triangulation_2.h>
  #include <CGAL/draw_triangulation_2.h>
  #include <CGAL/basic.h>
  #include <CGAL/Coercion_traits.h>
  #include <CGAL/IO/io.h>
  #include <fstream>
  typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
  typedef CGAL::Triangulation_2<K>                            Triangulation;
  typedef Triangulation::Point                                Point;

  template <typename A, typename B>
  typename CGAL::Coercion_traits<A,B>::Type
  binary_func(const A& a , const B& b){
      typedef CGAL::Coercion_traits<A,B> CT;
      CGAL_static_assertion((CT::Are_explicit_interoperable::value));
      typename CT::Cast cast;
      return cast(a)*cast(b);
  }

  int main(int argc, char**) {
    std::cout<< binary_func(double(3), int(5)) << std::endl;
    std::cout<< binary_func(int(3), double(5)) << std::endl;
    std::ifstream in("data/triangulation_prog1.cin");
    std::istream_iterator<Point> begin(in);
    std::istream_iterator<Point> end;
    Triangulation t;
    t.insert(begin, end);
    if(argc == 3) // do not test Qt5 at runtime
      CGAL::draw(t);
    return EXIT_SUCCESS;
   }
EOF
cat << EOF > "CMakeLists.txt"
  cmake_minimum_required(VERSION 3.1...3.15)
  find_package(CGAL COMPONENTS Qt5)
  add_definitions(-DCGAL_USE_BASIC_VIEWER -DQT_NO_KEYWORDS)
  add_executable(test main.cpp)
  target_link_libraries(test PUBLIC CGAL::CGAL_Qt5)
EOF
cd build
mytime cmake -DCGAL_DIR=$ROOT/config_dir ..
CGAL_PATH=$(cat CMakeCache.txt |grep -i cgal_dir)
if [ "$CGAL_PATH" != "CGAL_DIR:PATH=$ROOT/config_dir"]; then
  exit 1;
fi


make -j2
./test
cd ../..
rm -rf config_dir test_dir

#==-- CGAL_DIR pointing to the directory containing Git/CGALConfig.cmake (for the Git layout only) --==
mkdir test_dir
cd test_dir
touch main.cpp
mkdir build
cat << EOF > main.cpp
  #include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
  #include <CGAL/Triangulation_2.h>
  #include <CGAL/draw_triangulation_2.h>
  #include <CGAL/basic.h>
  #include <CGAL/Coercion_traits.h>
  #include <CGAL/IO/io.h>
  #include <fstream>
  typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
  typedef CGAL::Triangulation_2<K>                            Triangulation;
  typedef Triangulation::Point                                Point;

  template <typename A, typename B>
  typename CGAL::Coercion_traits<A,B>::Type
  binary_func(const A& a , const B& b){
      typedef CGAL::Coercion_traits<A,B> CT;
      CGAL_static_assertion((CT::Are_explicit_interoperable::value));
      typename CT::Cast cast;
      return cast(a)*cast(b);
  }

  int main(int argc, char**) {
    std::cout<< binary_func(double(3), int(5)) << std::endl;
    std::cout<< binary_func(int(3), double(5)) << std::endl;
    std::ifstream in("data/triangulation_prog1.cin");
    std::istream_iterator<Point> begin(in);
    std::istream_iterator<Point> end;
    Triangulation t;
    t.insert(begin, end);
    if(argc == 3) // do not test Qt5 at runtime
      CGAL::draw(t);
    return EXIT_SUCCESS;
   }
EOF
cat << EOF > "CMakeLists.txt"
  cmake_minimum_required(VERSION 3.1...3.15)
  find_package(CGAL COMPONENTS Qt5)
  add_definitions(-DCGAL_USE_BASIC_VIEWER -DQT_NO_KEYWORDS)
  add_executable(test main.cpp)
  target_link_libraries(test PUBLIC CGAL::CGAL_Qt5)
EOF
  cd build
  mytime cmake -DCGAL_DIR=$ROOT/CGALConfig.cmake ..
  if [ "$CGAL_PATH" != "CGAL_DIR:PATH=$ROOT"]; then
    exit 1;
  fi
  make -j2
  ./test
  cd ../..
  rm -rf config_dir test_dir
  #==-- configure all CGAL with -DWITH_examples=ON -DWITH_demos=ON -DWITH_tests=ON, and then launch CTest on a few labels. --==
  mkdir config_dir
  cd config_dir
  cmake -DWITH_examples=ON -DWITH_demos=ON -DWITH_tests=ON -DBUILD_TESTING=ON ..
  ctest -j2 -L AABB_tree
  cd ..
  rm -rf ./config_dir
    exit 0
  fi

  IFS=$old_IFS

  if [ -n "$TRAVIS_PULL_REQUEST_BRANCH" ] && [ "$ARG" != Polyhedron_demo ]; then
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
