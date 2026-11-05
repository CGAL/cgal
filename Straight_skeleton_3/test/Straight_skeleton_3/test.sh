#!/bin/bash

SCRIPT_NAME="${0##*/}"

# Print a concise summary of results (can be called from the trap)
print_summary() {
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
}

# When interrupted, print the summary and exit with 130 (standard for SIGINT)
print_summary_and_exit() {
  echo "Interrupted by user; stopping test script.";
  print_summary
  exit 130
}

trap 'print_summary_and_exit' INT TERM

BUILD_DIR="/home/mrouxell/git/CGAL/SLS3/Straight_skeleton_3/test/Straight_skeleton_3/build-release"
TIMEOUT_VALUE=120

SUCCESS=()
FAILURE=()
TIMEOUT=()
KILLED=()
IGNORED=()
POLYHEDRONS=$(\
  find /home/mrouxell/Data/Thingi10K/raw_meshes -type f -name '*.stl' | sort ;\
)

CURRENT_DATE=$(date +%Y-%m-%d)
CURRENT_TIME=$(date +%H:%M)
OUTPUT_DIRECTORY=test_results_${CURRENT_DATE}_${CURRENT_TIME}

cd $BUILD_DIR || exit 1

mkdir -p $OUTPUT_DIRECTORY

echo "Running tests with executable: $BUILD_DIR/test_skeleton_3"

for POLYHEDRON in $POLYHEDRONS; do
  CMD="./test_skeleton_3 $POLYHEDRON --save-times 1"

  echo ""
  echo "$CMD"
  # echo "$OUTPUT_FILE"

  # Ignore files above 1 MB
  # FILE_SIZE=$(stat -c%s "$POLYHEDRON")
  # if [ "$FILE_SIZE" -gt 1048576 ]; then
  #   echo "Skipping $POLYHEDRON (size: $FILE_SIZE bytes)"
  #   IGNORED+=("$CMD")
  #   continue
  # fi

  BASE_NAME=$(basename "$POLYHEDRON" | sed 's/\.[^.]*$//')
  OUTPUT_FILE=$OUTPUT_DIRECTORY/${BASE_NAME}.log

  timeout --foreground --signal=INT --kill-after=5s --preserve-status $TIMEOUT_VALUE $CMD > "$OUTPUT_FILE" 2>&1
  TIMEOUT_RES=$?
  echo "$TIMEOUT_RES"
  if [ "$TIMEOUT_RES" -eq 0 ]; then
    SUCCESS+=("$CMD")
  elif [ "$TIMEOUT_RES" -eq 143 ]; then
    TIMEOUT+=("$CMD")
  elif [ "$TIMEOUT_RES" -eq 137 ]; then
    KILLED+=("$CMD")
  else
    FAILURE+=("$CMD")
  fi
  echo ""
done

print_summary
