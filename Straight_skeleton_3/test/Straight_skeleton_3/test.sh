#!/bin/bash

SCRIPT_NAME="${0##*/}"

BUILD_DIR="/home/mrouxell/git/SLS3/Straight_skeleton_3/test/Straight_skeleton_3/build-debug"
TIMEOUT_VALUE=120

DATA_PATH=/home/mrouxell/git/SLS3/Straight_skeleton_3/test/Straight_skeleton_3/data

SUCCESS=()
FAILURE=()
TIMEOUT=()
POLYHEDRONS=$(\
  find $DATA_PATH/events-triangulated -type f -name '*.obj' | sort; \
  # find $DATA_PATH/events -type f -name '*.obj' | sort; \
  # find $DATA_PATH/mike-polytopes-triangulated -type f -name '*.obj' | sort ;\
  # find $DATA_PATH/wolfgang -type f -name '*.obj' | sort ;\
  # find $DATA_PATH/others -type f -name '*.obj' | sort; \
  # find $DATA_PATH/split-ev -type f -name '*.obj' | sort ;\
)

CURRENT_DATE=$(date +%Y-%m-%d)
CURRENT_TIME=$(date +%H:%M)
OUTPUT_DIRECTORY=test_results_${CURRENT_DATE}_${CURRENT_TIME}

cd $BUILD_DIR || exit 1

mkdir -p $OUTPUT_DIRECTORY

for POLYHEDRON in $POLYHEDRONS; do
  CMD="./test_skeleton_3 $POLYHEDRON"

  BASE_NAME=$(basename "$POLYHEDRON" | sed 's/\.[^.]*$//')
  OUTPUT_FILE=$OUTPUT_DIRECTORY/${BASE_NAME}.log

  echo ""
  echo "$CMD"
  echo "$OUTPUT_FILE"

  timeout --preserve-status $TIMEOUT_VALUE $CMD > "$OUTPUT_FILE" 2>&1
  TIMEOUT_RES=$?
  echo "$TIMEOUT_RES"
  if [ "$TIMEOUT_RES" -eq 0 ]; then
    SUCCESS+=("$CMD")
  elif [ "$TIMEOUT_RES" -eq 143 ]; then
    TIMEOUT+=("$CMD")
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
echo "TIMEOUT:"
for CMD in "${TIMEOUT[@]}"; do
  echo "$CMD"
done
