#!/bin/bash

if [ "$#" -lt 4 ]; then
    echo "Usage: $0 <input_file> <timeout> [component_params...]"
    exit 1
fi

timeout_bis() {
  timeout 5 sleep 10
}

INPUT_FILE=$1
TIMEOUT=$2
GRID_SIZE=$3
ERASE_ALL_DUPLICATE=$4

# Run with timeout, capture exit code
timeout "--foreground" "$TIMEOUT"s robustness_snap_polygon_soup "$INPUT_FILE" "$GRID_SIZE" "$ERASE_ALL_DUPLICATE"
EXIT_CODE=$?

# Interpret exit codes
declare -A TAGS
TAGS[0]="VALID_OUTPUT"
TAGS[1]="INPUT_IS_INVALID"
TAGS[2]="ROUNDING_FAILED"
TAGS[3]="SELF_INTERSECTING_OUTPUT"
TAGS[139]="SIGSEGV"
TAGS[11]="SIGSEGV"
TAGS[6]="SIGABRT"
TAGS[8]="SIGFPE"
TAGS[132]="SIGILL"
TAGS[124]="TIMEOUT"

TAG_NAME=${TAGS[$EXIT_CODE]:-UNKNOWN}
TAG_DESC=$([[ "$EXIT_CODE" -eq 0 ]] && echo "OK" || echo "Error")

# Output JSON
echo "{\"TAG_NAME\": \"$TAG_NAME\", \"TAG\": \"$TAG_DESC\"}"