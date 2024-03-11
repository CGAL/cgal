#!/bin/bash

BUILD_DIR="build"


test_2d () {
  cd "$BUILD_DIR" || exit 1
  SUCCESS=()
  FAILURE=()
  for POLYGON_ID in $(seq 12); do
    CMD="./StraightSkel 2d $POLYGON_ID --no-window"
    echo ""
    echo "$CMD"
    $CMD
    if [ $? -eq 0 ]; then
      SUCCESS+=("$CMD")
    else
      FAILURE+=("$CMD")
    fi
    echo ""
  done
  echo ""
  echo "SUCCESS:"
  for CMD in "${SUCCESS[@]}"; do
    echo "$CMD"
  done
  echo ""
  echo "FAILURE:"
  for CMD in "${FAILURE[@]}"; do
    echo "$CMD"
  done
  echo ""
}


test_3d () {
  cd "$BUILD_DIR" || exit 1
  SUCCESS=()
  FAILURE=()
  POLYHEDRONS=$(find ../res/polyhedrons/events -type f -name '*.obj' | sort)
  for POLYHEDRON in $POLYHEDRONS; do
    CMD="./StraightSkel 3d load $POLYHEDRON --no-window"
    echo ""
    echo "$CMD"
    $CMD
    if [ $? -eq 0 ]; then
      SUCCESS+=("$CMD")
    else
      FAILURE+=("$CMD")
    fi
    echo ""
  done
  echo ""
  echo "SUCCESS:"
  for CMD in "${SUCCESS[@]}"; do
    echo "$CMD"
  done
  echo ""
  echo "FAILURE:"
  for CMD in "${FAILURE[@]}"; do
    echo "$CMD"
  done
  echo ""
}


case "$1" in
'2d')
  test_2d
  ;;
'3d')
  test_3d
  ;;
*)
  test_3d
esac
