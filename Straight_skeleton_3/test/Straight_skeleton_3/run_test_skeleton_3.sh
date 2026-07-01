#!/bin/bash

GREEN="\\033[1;32m"
NORMAL="\\033[0;39m"
RED="\\033[1;31m"
TIMEOUT_VALUE=1800

cleanup() {
  echo -e "\nInterrupted by user"
  pkill -P $$ 2>/dev/null
  exit 1
}
trap cleanup SIGINT SIGTERM

line_number=0
while IFS= read -r line; do
  ((line_number++))
  [ -z "$line" ] && continue

  expanded_line=$(eval echo "$line")
  f=$(echo "$expanded_line" | awk '{print $1}')
  echo "[$line_number] ==== $line"

  # Skip if expanded line is empty
  [ -z "$expanded_line" ] && continue

  echo "[$line_number] ==== $expanded_line"

  CMD="build-release/test_skeleton_3 $expanded_line"

  if timeout --preserve-status --foreground $TIMEOUT_VALUE $CMD; then
    echo -e " ==> $GREEN SUCCEED"
    echo -e -n "$NORMAL"
  else
    STATUS=$?
    if [ $STATUS -eq 124 ]; then
      echo -e " ==> $RED TIMEOUT"
    else
      echo -e " ==> $RED FAILED (status: $STATUS)"
    fi
    echo -e -n "$NORMAL"
  fi
done < test_skeleton_3.cmd