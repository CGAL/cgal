#!/bin/sh
# out-of-source build

BUILD_DIR="build"
DB_2D="skeldata2d.db3"

# exit on error
set -e

build_dir () {
  mkdir -p "$BUILD_DIR"
  cd "$BUILD_DIR" || exit 1
}

build () {
  cmake -D CMAKE_BUILD_TYPE=Debug ..
  make
}

build_db_2d () {
  ./StraightSkel 2d 0 > /dev/null || true   # Create initial database
  cat ../res/create_polygons.sql | sqlite3 "$DB_2D"
}

clean () {
  rm -rf "$BUILD_DIR"
  rm -rf CMakeFiles
  rm -rf CMakeCache.txt
  rm -rf *.cmake
}


case "$1" in
'clean')
  clean
  ;;
'db')
  build_dir
  rm -rf "$DB_2D"
  build_db_2d
  ;;
*)
  build_dir
  build
  if [ ! -f "$DB_2D" ]; then
    build_db_2d
  fi
  ;;
esac
